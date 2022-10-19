# this python file is for using rdkit generated basic structure info to generate descriptors
from get_Info import *






#INPUT_FILES_______________________________________________________________________________________________________

EN_data_path = "__data_empirical/EN.txt"
with open(EN_data_path,"r") as input:
    EN_dict = {}
    count_ = 0
    for row in input.read().split("\n"):
        if count_ == 0:
            count_ += 1
            continue
        elem , value = row.split(" ")
        EN_dict[elem] = float(value)\






#FUNCTIONS_______________________________________________________________________________________________________

# find the index of donor atom automatically
def get_LigAcType2LigAcId(lig_ac_type_str,smiles):
    # check whether donor atom is carbene or not
    locate_carbene = True if lig_ac_type_str == "C" else False

    # change str to list of individual atom str
    try:
        lig_ac_type_list = lig_ac_type_str.split(";")
    except:
        print("smiles:",smiles," have no donor type, pls add manually in excel")
        return None

    mol = Chem.MolFromSmiles(smiles)
    lig_ac_id_list = []
    for atom in mol.GetAtoms():
        a_id = atom.GetIdx()
        elem = atom.GetSymbol()
        # check if the element match donor atom
        if elem in lig_ac_type_list:
            # if there is carbene, check for 2 neighbouring N atoms to confirm
            if locate_carbene:
                nghbr_N_count = 0
                for atom in atom.GetNeighbors():
                    if atom.GetSymbol() == "N":
                        nghbr_N_count+=1
                # skip and search for the next atom if cannot find carbene
                if nghbr_N_count != 2:
                    continue
            # store a_id of donor atom
            lig_ac_id_list.append(a_id)
    # return error msg if one more match exists
    if len(lig_ac_id_list) != len(lig_ac_type_list):
        msg = "ERROR: there maybe more than 1 matches for donor atom, please locate donor atoms manually - \n  lig_ac_type_str= "+lig_ac_type_str+"\n  smiles= "+smiles
        return msg
    # return sorted target id list
    else:
        return sorted(lig_ac_id_list)


# get the information about donor atom of the ligand
def get_Info2LigCaInfo(info_dict, output_format = None):
    lig_ac_info_dict = {}
    for current_target, value in info_dict.items():
        # get target_atom information
        target_atom_info = value['target_atom_info']

        # get electronegativity
        EN = EN_dict[target_atom_info['element']]

        # get hybridization
        if target_atom_info['hybrid'] == 'S':
            hybrid = 1
        elif target_atom_info['hybrid'] == 'SP':
            hybrid = 0.5
        elif target_atom_info['hybrid'] == 'SP2':
            hybrid = 0.33
        elif target_atom_info['hybrid'] == 'SP3':
            hybrid = 0.25
        else:
            hybrid = 0

        # get aromatic ring size
        arom = 0
        for ring in target_atom_info['ring']:
            if ring['ring_aromaticity']:
                arom = ring['ring_size']

        # integrate the information together
        lig_ac_info_dict[current_target] = { 
            'EN': EN,
            'hybrid': hybrid,
            'aromaticity':  arom,} # arom = invovled Ar bond count

    # change output format to list form accordin to order of input info_dict
    if output_format == "list":
        lig_ac_info_list = []
        for lig_ac in lig_ac_info_dict:
            lig_ac_info_list.append(list(lig_ac_info_dict[lig_ac].values()))
        return lig_ac_info_list
    
    # return ligand ac info dict
    return lig_ac_info_dict


# get the sigma donating descriptor
def get_Info2SigmaFeat(info_dict, cutoffRadius = None):
    sigma_score_total = 0
    for current_target,value in info_dict.items():
        sigma_score = 0

        # get the target atom info
        target_atom_info = value['target_atom_info']
        # get info of every other atoms
        for atom_info in value['atom_info_list']:
            
            # check the EN and cal difference of target atom to current atom
            EN_diff = EN_dict[target_atom_info['element']] - EN_dict[atom_info['element']]
            
            # normalize the difference !!! Normalization
            EN_diff = EN_diff/4

            # check hybridization state of atom
            if atom_info['hybrid'] == 'S':
                hybrid = 1
            elif atom_info['hybrid'] == 'SP':
                hybrid = 0.5
            elif atom_info['hybrid'] == 'SP2':
                hybrid = 0.33
            elif atom_info['hybrid'] == 'SP3':
                hybrid = 0.25
            else:
                hybrid = 0

            # get bond displ from target atom
            bond_displ = atom_info['bonddisplacement']

            # if there is a cut off radius (displ), all atom beyond the displ will not be included
            if cutoffRadius != None:
                if int(bond_displ) > cutoffRadius:
                    continue

            # +ve means electron donating; -ve means electron withdrawing
            sigma_score +=  EN_diff * hybrid / bond_displ**2
            #print(sigma_score)
        sigma_score_total += sigma_score

    # cal the average of score be dividing per target atoms
    sigma_score_avg = sigma_score_total/len(info_dict.keys())
    return sigma_score_avg
    

# total steric: reference Topological steric effect index and its application
def get_Info2SterWholeFeat(info_dict):
    whole_steric_score_total = 0
    for current_target, value in info_dict.items():
        TSEI = 0 # topological distance 

        # get the target atom info
        target_atom_info = value['target_atom_info']
        # get info of every other atoms
        for atom_info in value['atom_info_list']:
            # check the bond displacement
            bond_displ = atom_info['bonddisplacement']
            # cal TSEI for current atom and add to total
            TSEI += 1/ bond_displ**3
    
        # total number of bond
        total_bond_count = value["Bond_Overall_info"]["total_bond_count"]
        Ar_bond_count = value["Bond_Overall_info"]["Ar_bond_count"]
        NonArRing_bond_count = value["Bond_Overall_info"]["NonArRing_bond_count"]

        if total_bond_count != 0:
            # cal the bond count as percentage
            rot_able_bond_percent = (total_bond_count - Ar_bond_count - 0.5*NonArRing_bond_count)/total_bond_count
        else:
            rot_able_bond_percent = 0
        
        # use TSEI times the average of (rot_able_bond_percent + 1) so that unrotatble bond won't be zero
        whole_steric_score = TSEI * (rot_able_bond_percent + 1)/2

        # add to the total steric score
        whole_steric_score_total += whole_steric_score

    # cal the average of score by dividing per target atoms
    whole_steric_score_avg = whole_steric_score_total / len(info_dict.keys())
    return whole_steric_score_avg



#EXAMPLES_FOR_CHECKING____________________________________________________________________________________________
