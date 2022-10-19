import pandas as pd
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw





#INPUT_FILES_______________________________________________________________________________________________________

default_output_path = '_test/mol.o.png'






#FUNCTIONS_______________________________________________________________________________________________________

# to allocate task for data type
def get_AllInfo(input_,target_ids = None,IsFragment = None,data_type = "smiles"):
    if data_type == "smiles":
        mol = Chem.MolFromSmiles(input_)
        #show_Mol2AutoAmapImg(mol)
        return get_Mol2AllInfo(mol,target_ids=target_ids,IsFragment=IsFragment)
    # if data_type is a json
    elif data_type == "json":
        info = {}
        for key,value in json.loads(input_).items():
            # if the input is a list
            if type(value) == list:
                for smil in value:
                    mol = Chem.MolFromSmiles(smil)
                    mol_info = get_Mol2AllInfo(mol,target_ids=target_ids,IsFragment=IsFragment)
                info[key] = mol_info
            # if the input is a list
            else:
                smil = value
                mol = Chem.MolFromSmiles(smil)
                #show_Mol2AutoAmapImg(mol)
                mol_info = get_Mol2AllInfo(mol,target_ids=target_ids,IsFragment=IsFragment)
                info[key] = mol_info
        # return info as a dict
        return info


# show a automatically atom-mapped graph from mol
def show_Mol2AutoAmapImg(mol,output_path = default_output_path):
    print("Please look at image at ",output_path,"and input_ target_ids before proceeding e.g. 1 or 1,2,3:")
    atom_info_list = []
    for atom in mol.GetAtoms():
        a_id = atom.GetIdx()
        elem = atom.GetSymbol()
        atom.SetAtomMapNum(a_id)
    Draw.MolToFile(mol,output_path)    


# !!!! becareful, only works with mol generated from smiles
# changed reference from https://www.rdkit.org/docs/Cookbook.html#count-ring-systems
def get_Mol2RingInfo(mol, includeSpiro=False, Mapped = False,aid2amap_dict = {}):
    # get a set of aromatic atoms
    Ar_atom_set = set([i.GetIdx() for i in mol.GetAromaticAtoms()])
    # get rings from rdkit
    ri = mol.GetRingInfo()

    systems = []
    for ring in ri.AtomRings():
        # all a_id in that ring
        ringats = set(ring)
        
        # check if ring is aromatic
        arom = True if ringats == Ar_atom_set.intersection(ringats) else False

        # copied code from reference
        nSystems = []
        for system in systems:
            nInCommon = len(ringats.intersection(system['ringAts']))
            if nInCommon and (includeSpiro or nInCommon>1) and (system['aromaticity'] == arom):
                ringats = ringats.union(system['ringAts'])
            else:
                nSystems.append(system)

        # give ring a_ids and arom info
        nSystems.append({'ringAts':ringats,'aromaticity':arom,'ringSize':len(ringats)})
        systems = nSystems
    
    # changing name only
    ring_info_list = systems

    # change a_id to amap and add ringsize (cannot be added directly to above code as it depends on identifying the a_id)
    if Mapped:
        for ring_info in ring_info_list:
            list_ = list(ring_info['ringAts'])
            for i in range(len(list_)):
                list_[i] = aid2amap_dict[list_[i]]
            ring_info['ringAts'] = list(set(list_))
    return ring_info_list


# get a full info from mol
def get_Mol2AllInfo(mol,target_ids = None,IsFragment = False,Mapped = False): # target_ids must be interger
    # !!!! note that target_ids = atom map if Mapped = True; 
    #                target_ids = atom id if Mapped = False
    if target_ids == None and IsFragment == None:
        print("ERROR: please specify the task")
        return

    # turn to list for convenience 
    target_list = target_ids
    if type(target_ids) != list:
        target_list = [target_ids]

    # only used when IsFragment == True
    brkpt = None

    aid2amap_dict = {}
    amap2aid_dict = {}
    if IsFragment or Mapped:
        for atom in mol.GetAtoms():
            a_id = atom.GetIdx()
            amap = atom.GetAtomMapNum()
            elem = atom.GetSymbol()
            # generate the a dict for a_id and a_map for convience 
            if Mapped:
                aid2amap_dict[a_id] = amap
                amap2aid_dict[amap] = a_id
            # locate the fragment position as break pt
            if IsFragment and elem == "*": 
                brkpt = atom
                target_list = [atom.GetIdx()]
    
    output_dict = {}
    target_atom_info = {}
    for current_target in target_list:
        bond_displ_matrix = AllChem.GetDistanceMatrix(mol)
        ring_info_list = get_Mol2RingInfo(mol)
        # change ring_info_list to amap instead of aid
        if Mapped:
            for ring_info in ring_info_list:
                list_ = list(ring_info['ringAts'])
                for i in range(len(list_)):
                    list_[i] = aid2amap_dict[list_[i]]
                ring_info['ringAts'] = set(list_)

        atom_info_list = []
        for atom in mol.GetAtoms():
            a_id = atom.GetIdx()
            # If not mapped, assume atom mapping is the same as atom id
            if not Mapped:
                atom.SetAtomMapNum(a_id)

            # property of atom
            elem = atom.GetSymbol()
            arom =  atom.GetIsAromatic()
            hybrid = atom.GetHybridization()
            radicale = atom.GetNumRadicalElectrons()

            # get bond displacement 
            current_target_id = amap2aid_dict[current_target] if Mapped else current_target
            bond_displ = bond_displ_matrix[current_target_id][a_id]

            # generate output for each atom
            atom_info = {
                        'aid':  a_id,
                        'amap': aid2amap_dict[a_id] if Mapped else -1,
                        'element': elem,
                        'aromaticity': arom,
                        'hybrid': str(hybrid),
                        'radical_e': radicale,
                        'bonddisplacement': bond_displ
                        }

            # target atom info are stored separately from others
            if a_id == current_target_id: 
                # generate ring info for target atom
                ring = []
                for ring_info in ring_info_list:
                    if current_target in ring_info['ringAts']:
                        ring.append({   'ring_size':ring_info['ringSize'],
                                        'ring_aromaticity':ring_info['aromaticity']})
                atom_info['ring'] = ring

                target_atom_info = atom_info
                continue

            # push output info for each atom except target atoms
            atom_info_list.append(atom_info)

        # generate bond info
        total_bond_count = 0
        Ar_bond_count = 0
        NonArRing_bond_count = 0
        for bond in mol.GetBonds():
            if brkpt != None:
                # skip counting the break pt bond if isFragment (ignore the *-C bond)
                if brkpt.GetIdx() in [bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()]:
                    continue
            total_bond_count+=1
            if bond.GetIsAromatic():
                Ar_bond_count+=1
            elif bond.IsInRing():
                NonArRing_bond_count+=1
        Bond_Overall_info = {
                            'total_bond_count': total_bond_count,
                            'Ar_bond_count': Ar_bond_count,
                            'NonArRing_bond_count': NonArRing_bond_count
                            }

        # finally, add all the info together
        final_info = {
                            "target_atom_info":target_atom_info,
                            'atom_info_list':atom_info_list,
                            'Bond_Overall_info' : Bond_Overall_info,
                            'Ring_info_list':ring_info_list
                            }
        
        if IsFragment:
            output_dict = final_info
            break
        # return a dict with current_target if it is not fragment
        #CHANGETO if current_target is a single input_ return the final_info directly instead
        else:
            output_dict[current_target] = final_info
    return output_dict



# get only neighbor info from mol
def get_Mol2NghbrInfo(mol,bond_atom_exception_list = [],Mapped = False,target_ids = None,IsFragment = False): # target_ids must be interger
    if target_ids == None and IsFragment == None:
        print("ERROR: please specify the task")
        return

    # turn to list for convenience 
    target_list = target_ids
    if type(target_ids) != list:
        target_list = [target_ids]

    brkpt = None #only used when IsFragment == True
    aid2amap_dict = {}
    amap2aid_dict = {}
    if IsFragment or Mapped:
        for atom in mol.GetAtoms():
            a_id = atom.GetIdx()
            amap = atom.GetAtomMapNum()
            elem = atom.GetSymbol()
            if Mapped:
                aid2amap_dict[a_id] = amap
                amap2aid_dict[amap] = a_id
            if IsFragment and elem == "*":  #the fragment position
                brkpt = atom
                target_list = [atom.GetIdx()]
    
    output_dict = {}
    target_atom_info = {}
    for current_target in target_list:
        bond_displ_matrix = AllChem.GetDistanceMatrix(mol)
        ring_info_list = get_Mol2RingInfo(mol)
        # change ring_info_list to mapped form
        if Mapped:
            for ring_info in ring_info_list:
                list_ = list(ring_info['ringAts'])
                for i in range(len(list_)):
                    list_[i] = aid2amap_dict[list_[i]]
                ring_info['ringAts'] = set(list_)

        #above_is_same_as_get_Mol2AllInfo___________________________________________________________________________________
        current_atom = mol.GetAtomWithIdx(amap2aid_dict[current_target]) if Mapped else mol.GetAtomWithIdx(a_id)
        nghbr_atom_info_list = []
        # get info of neighbors atom as well as the target atom
        for atom in current_atom.GetNeighbors()+(current_atom,):
            a_id = atom.GetIdx()
            
            # this id is changed to amap form if Mapped
            current_nghbr_id = aid2amap_dict[a_id] if Mapped else a_id
            # skip if in exception_list
            if current_nghbr_id in bond_atom_exception_list:
                continue

            # If not mapped, assume atom mapping is the same as atom id
            if not Mapped:
                atom.SetAtomMapNum(a_id)

            # property of neighbor atom
            elem = atom.GetSymbol()
            arom =  atom.GetIsAromatic()
            hybrid = atom.GetHybridization()
            radicale = atom.GetNumRadicalElectrons()

            # get bond info of neighbor atom
            bonded_atom_info_list = []
            bond_class_str_list = []
            for bond in atom.GetBonds():
                bonded_atom = bond.GetOtherAtom(atom)

                # not include the target_atom
                if bonded_atom.GetIdx() == current_atom.GetIdx():
                    continue
                bond_type = bond.GetBondTypeAsDouble()

                # property of bonded atom of neighbor atom
                bonded_atom_id = bonded_atom.GetIdx()
                bonded_atom_elem = bonded_atom.GetSymbol()
                bonded_atom_arom =  bonded_atom.GetIsAromatic()
                bonded_atom_hybd = bonded_atom.GetHybridization()
                bonded_atom_info = {
                        'aid': bonded_atom_id,
                        'amap': aid2amap_dict[bonded_atom_id] if Mapped else -1,
                        'element': bonded_atom_elem,
                        'aromaticity': bonded_atom_arom,
                        'hybrid': str(bonded_atom_hybd),
                        'radical_e': radicale,
                        'bond_type': bond_type,
                        }
                bonded_atom_info_list.append(bonded_atom_info)

                if bond_type == 1.0:
                    bond_text = "-"
                elif bond_type == 1.5:
                    bond_text = ":"
                elif bond_type == 2.0:
                    bond_text = "="
                else:
                    bond_text = "?"
                if bonded_atom_arom:
                    bond_Ar_text = "Ar"
                else:
                    bond_Ar_text = "NonAr"

                # generate a bond class str for the bonded atom - neighbour atom pair
                bond_class_str = elem + bond_text + bonded_atom_elem + "("+bond_Ar_text+";"+str(bonded_atom_hybd)+")"
                bond_class_str_list.append(bond_class_str)

            # integrate the above info in atom_info for a neighbor atom
            atom_info = {
                        'aid': a_id,
                        'amap': aid2amap_dict[a_id] if Mapped else -1,
                        'element': elem,
                        'aromaticity': arom,
                        'hybrid': str(hybrid),
                        'radical_e': radicale,
                        'bond_list':bond_class_str_list,
                        'bonded_atom_info_list':bonded_atom_info_list,
                        }


             # target atom info saved separately
            if current_nghbr_id == current_id:  
                ring = []
                for ring_info in ring_info_list:
                    if current_target in ring_info['ringAts']:
                        ring_size = len(ring_info['ringAts'])
                        ring_arom = ring_info['aromaticity']
                        ring.append({'ring_size':ring_size,'ring_aromaticity':ring_arom})
                atom_info['ring'] = ring
                target_atom_info = atom_info
                continue
            nghbr_atom_info_list.append(atom_info)

        # integrate all info
        final_info = {
                            "target_atom_info":target_atom_info,
                            'neighbor_atom_info_list':nghbr_atom_info_list,
                            'Ring_info_list':ring_info_list
                            }
        
        if IsFragment:
            output_dict = final_info
            break
        # return a dict with current_target if it is not fragment
        #CHANGETO if current_target is a single input_ return the final_info directly instead
        else:
            output_dict[current_target] = final_info
    return output_dict



#EXAMPLES_FOR_CHECKING____________________________________________________________________________________________
