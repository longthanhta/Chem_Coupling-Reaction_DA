import pandas as pd
from rdkit.Chem import AllChem as Chem
'''from rdkit import Chem
from rdkit.Chem import AllChem'''
from rdkit.Chem import Draw
import re
import rdkit.Chem.rdFMCS as MCS


#INPUT_FILES_______________________________________________________________________________________________________


#FUNCTIONS_______________________________________________________________________________________________________

# create an img file that has rxn smarts
#CHANGETO create_Smat2Img
def Draw_temp(smat,filename,skip_error=False):
    #CHANGETO Do_create_Smat2Img
    def Do_Draw_temp(smat,filename):
        # output an img from rxn smarts
        rxn = Chem.ReactionFromSmarts(smat)
        rxn_img = Draw.ReactionToImage(rxn)
        rxn_img.save(filename)

    # ignore error or not
    if skip_error:
        try:
            Do_Draw_temp(smat,filename)
        except:
            pass
    else:
        Do_Draw_temp(smat,filename)

# get bond type symbol from bond
#CHANGETO get_Bond2BondType
def get_bond_type(bond):
    bond_type=str(bond.GetBondType())
    if bond_type=='SINGLE':
        bond_type='-'
    if bond_type=='DOUBLE':
        bond_type='='
    if bond_type=='TRIPLE':
        bond_type='#'
    if bond_type=='AROMATIC':
        bond_type=':'
    return bond_type

# get hybrd dict based on amap
#CHANGETO get_Smil2Hybd
def get_hybridization(smil):
    # when elem == C or N
    elem_to_get = ['C','N']
    mol = Chem.MolFromSmiles(smil)
    hybd = ['' for atom in mol.GetAtoms()]

    for atom in mol.GetAtoms():
        # does not care if elem not = C or N
        if atom.GetSymbol() not in elem_to_get: continue
        adj_bonds = []
        a_id = atom.GetIdx()

        # record bond types
        for bond in mol.GetBonds():
            a_id1 = bond.GetBeginAtomIdx()
            a_id2 = bond.GetEndAtomIdx()
            if a_id == a_id1 or a_id == a_id2:
                adj_bonds.append(get_bond_type(bond))

        # record bond types and corresponding count
        result = pd.value_counts(adj_bonds)
        # assign hybrd
        if ':' in result.index:
            hybd[a_id] = 'sp2'
        elif '#' in result.index:
            hybd[a_id] = 'sp1'
        elif '=' in result.index:
            if result['='] == 2:
                hybd[a_id] = 'sp1'
            else:
                hybd[a_id] = 'sp2'
        else:
            hybd[a_id] = 'sp3'
    return hybd

# get a bond class str from input
#CHANGETO get_BondClassStr
def get_bond_str(hybd_atom1,hybd_atom2,a_map1,a_map2,elem1,elem2,no_H1,no_H2,bond_type):
    atom_str1=elem1+hybd_atom1+'H'+str(no_H1)+'_'+str(a_map1)
    atom_str2=elem2+hybd_atom2+'H'+str(no_H2)+'_'+str(a_map2)
    if elem1+str(a_map1)>elem2+str(a_map2):
        bond_str=atom_str1+bond_type+atom_str2
    else:
        bond_str=atom_str2+bond_type+atom_str1
    return bond_str

# get the a_ids of Reaction Center and Leaving Group Center atom
#CHANGETO get_Smil2RcAId_AND_LvgCaAId
def get_AM_cut_index(smil):
    mol = Chem.MolFromSmiles(smil)
    rc_id_list=[]
    lvg_ca_id_list=[]
    for bond in mol.GetBonds():
        # properties
        a_id1 = bond.GetBeginAtomIdx()
        a_id2 = bond.GetEndAtomIdx()

        atom1=mol.GetAtomWithIdx(a_id1)
        atom2=mol.GetAtomWithIdx(a_id2)

        amap1=atom1.GetAtomMapNum()
        amap2=atom2.GetAtomMapNum()
        # if amap = 0 -> lvg group atom
        if amap1==0 and amap2!=0:
            rc_id_list.append(a_id2)
            lvg_ca_id_list.append(a_id1)
        elif amap2==0 and amap1!=0:
            rc_id_list.append(a_id1)
            lvg_ca_id_list.append(a_id2)
    return list(set(rc_id_list)), list(set(lvg_ca_id_list))

# get list of bond class from smiles
# !!! do not treat unmapped mol
#CHANGETO get_Smil2BondClassList
def get_list_of_bond(smil):
    # bond has at least 1 atom map
    bond_class_list=[]
    # Check unmapped molecule, if there is a molecule without any mapping atom, ignore it
    unmapped_molecule=True
    mol = Chem.MolFromSmiles(smil)
    hybrid = get_hybridization(smil)

    for bond in mol.GetBonds():
        #bond properties
        a_id1 = bond.GetBeginAtomIdx()
        a_id2 = bond.GetEndAtomIdx()

        bond_type=get_bond_type(bond)

        atom1=mol.GetAtomWithIdx(a_id1)
        atom2=mol.GetAtomWithIdx(a_id2)

        no_H1=atom1.GetTotalNumHs()
        no_H2=atom2.GetTotalNumHs()

        amap1=atom1.GetAtomMapNum()
        amap2=atom2.GetAtomMapNum()

        hybrid_atom1 = hybrid[a_id1]
        hybrid_atom2 = hybrid[a_id2]

        elem1=atom1.GetSymbol()
        elem2=atom2.GetSymbol()
        # change halogen symbol to X
        if elem1 in ['Cl','Br','F','I']:
            elem1='X'
        if elem2 in ['Cl','Br','F','I']:
            elem2='X'
        
        bond_class_str=get_bond_str(hybrid_atom1,hybrid_atom2,amap1,amap2,elem1,elem2,no_H1,no_H2,bond_type)
        if not (amap1 ==0 and amap2==0):
            unmapped_molecule=False
            bond_class_list.append(bond_class_str)
    if unmapped_molecule:
        return []
    return bond_class_list

# get bond that is different from rxn
#CHANGETO get_rxn2DiffBond
def find_diff_bond(rxn):

    rct,prd = rxn.split('>>')
    # get list of bonds class
    rct_bond_class_list=list(set((get_list_of_bond(rct))))
    prd_bond_class_list=list(set((get_list_of_bond(prd))))

    rct_broken_list = []
    rct_changed_list = []
    prd_formed_list = []
    prd_changed_list = []
    rcts = rct.split('.')
    prds = prd.split('.')

    # get bond broken and bond changed for rct and prd
    for rct in rcts:
        bond_class_list = get_list_of_bond(rct)
        broken_list,changed_list = bond_subtraction(bond_class_list, prd_bond_class_list)
        rct_broken_list.append(broken_list)
        rct_changed_list.append(changed_list)
    for prd in prds:
        bond_class_list = get_list_of_bond(prd)
        formed_list, changed_list = bond_subtraction(bond_class_list, rct_bond_class_list)
        prd_formed_list.append(formed_list)
        prd_changed_list.append(changed_list)

    return rct_broken_list, rct_changed_list, prd_formed_list, prd_changed_list

# remove atom map from bond class list
#CHANGETO remove_Amap5BondClassList
def remove_AM_list(bond_class_list):
    bond_class_list_output=[]
    for bond in bond_class_list:
        bond_class_list_output.append(re.sub(r"_\d+","",bond))
    return bond_class_list_output

# remove atom map from bond class str
#CHANGETO remove_Amap5BondClassStr
def remove_AM_single(bond_class_str):
    return re.sub(r"_\d+","",bond_class_str)

# sort bond class str with atom map
#CHANGETO sort_BondClassStr
def sort_bond(bond_class_str):
    # identify bond type
    for bond_type in ['-','=','#',':']:
        if bond_type in bond_class_str:
            break
    # locate bond type
    bond_type_pos = bond_class_str.index(bond_type)
    # get string of atom1 and atom2
    atom1_str = bond_class_str[:bond_type_pos]
    atom2_str = bond_class_str[bond_type_pos+1:]
    # sort string of atom1 and atom2
    [atom1_str,atom2_str] = sorted([atom1_str,atom2_str])
    return atom1_str + bond_type + atom2_str

# sort bond class without atom map with alphabet order
#CHANGETO sort_BondClassWOAmap
# accept 2 types of data: 1. list ;2. list_of_list
def sort_bond_without_AM(input_list,list_of_list=True,changing_bond=False):
    # when input_list is a simple list
    if not list_of_list:
        output_sorted=[]
        for bond_class_str in input_list:
            # sort a single bond class string, remove amap part in process
            output_sorted.append(sort_bond(re.sub(r"_\d+","",bond_class_str)))
        # sort the whole output list
        output_sorted.sort()
        return output_sorted

    # when input_list is a list of list
    if changing_bond:
        output_sorted=[]
        for bond_prdair in input_list:
            bond_list=bond_prdair.split('>')
            # sort a changing bond class string, remove amap part in process
            output_sorted.append(sort_bond(re.sub(r"_\d+","",bond_list[0]))+'>'+sort_bond(re.sub(r"_\d+","",bond_list[1])))
        # sort the whole output list
        output_sorted.sort()
        return output_sorted   
    else:
        output_sorted=[]
        # for each molecule in input_list
        for mol in input_list:
            # for each molecule in input_list
            mol_bond_class_list=[]
            for bond_class_str in mol:
                # sort a changing bond class string, remove amap part in process
                mol_bond_class_list.append(sort_bond(re.sub(r"_\d+","",bond_class_str)))
            # sort the output within mol level
            mol_bond_class_list.sort()
            output_sorted.append(mol_bond_class_list)
        # sort the whole output list
        output_sorted.sort()
        return output_sorted       

# remove bond type from bond class string
#CHANGETO remove_BondType5BondClassStr
def remove_bond_type(bond_class_str):
    # identify bond type
    for bond_type in ['-','=','#',':']:
        if bond_type in bond_class_str:
            break
    # return the 2 atoms in the bond class string
    atom_str_list = bond_class_str.split(bond_type)
    # modify the atom string
    atom_str_list = [atom.split('_')[0].split('H')[0].split('sp')[0]+'_'+atom.split('_')[1] for atom in atom_str_list]
    # sort the atom string list
    atom_str_list.sort()
    # form back bond class string without bond type
    bond_class_wo_bond_type_str = '-'.join(atom_str_list)
    return bond_class_wo_bond_type_str

# find the difference of bond between two bond list
#CHANGETO subtract_Bond
def bond_subtraction(bond_class_list1, bond_class_list2,debug_mode=False): # Set False for debug mode
    # get difference of bond as a list
    diff_bond_class_list = list(set(bond_class_list1) - set(bond_class_list2))
    # get remove bond type version of bond_class_list2 for formalising the bond type to faciliate later matching
    bond_class_wo_bond_type_list2 = [remove_bond_type(bond_class_str) for bond_class_str in bond_class_list2]
    
    new_bond_list = []
    changed_bond_list = []
    for bond_class_str in diff_bond_class_list:
        # get changed bond
        if remove_bond_type(bond_class_str) in bond_class_wo_bond_type_list2:
            changed_bond_list.append(bond_class_str)
        # get formed bond
        else:
            new_bond_list.append(bond_class_str)
    if debug_mode:
    		new_bond_list = remove_AM_list(new_bond_list)
    return new_bond_list, changed_bond_list

# get the total No. of atoms that do not have atom map
#CHANGETO smil2NonAmapTotal
def get_noAM_number(smil):
    cnt = 0
    mol = Chem.MolFromSmiles(smil)
    for atom in mol.GetAtoms():
        # if no atom map, count +1
        if atom.GetAtomMapNum() == 0:
            cnt += 1 
    return cnt

# check whether 2 smiles are the same or not structurally
#CHANGETO checkidentical_Smil
def check_identical_Smiles(smil1, smil2):
    mol1 = Chem.MolFromSmiles(smil1)
    mol2 = Chem.MolFromSmiles(smil2)
    try:
        # use rdkit HasSubstructMatch to find whether they are both within each other's structure
        return mol2.HasSubstructMatch(mol1) and mol1.HasSubstructMatch(mol2)
    except:
        return False

# get the smiles of non atom mapped group
#CHANGETO get_SmilList2SmilListWOAmap
def get_no_AM_group(smil_list):
    smil_list_out = []
    for current_smil in smil_list:
        # get the reaction center atom id list and lvg group center atom id list
        rc_id_list, lvg_ca_id_list = get_AM_cut_index(current_smil)
        smil_list_wo_amap=[]
        smil_list_w_amap=[]
        
        # cut out the scfd
        for lvg_ca in lvg_ca_id_list:
            # get a editable mol
            rwmol=Chem.RWMol(Chem.MolFromSmiles(current_smil))
            # remove lvg group center atom by atom id
            rwmol.RemoveAtom(lvg_ca)
            # convert rwmol to smiles list
            smil_list=Chem.MolToSmiles(rwmol).split('.')
            for smil in smil_list:
                # check if the smiles contains amap by using ':No.'
                contain_amap=False
                if re.search(':\d+', smil):
                    contain_amap=True
                # identify if it is a completely separate scfd
                if contain_amap:
                    smil_list_w_amap.append(smil)
        # cut out leaving group 
        for rc in rc_id_list:
            # get a editable mol
            rwmol=Chem.RWMol(Chem.MolFromSmiles(current_smil))
            # remove lvg reaction center atom by atom id
            rwmol.RemoveAtom(rc)    
            # convert rwmol to smiles list
            smil_list=Chem.MolToSmiles(rwmol).split('.')
            for smil in smil_list:
                contain_amap=False
                # check if the smiles contains amap by using ':No.'
                if re.search(':\d+', smil):
                    contain_amap=True
                # identify if it is a completely separate leaving group
                if not contain_amap:
                    smil_list_wo_amap.append(smil)

        # looping scfd and leaving group to check identical or not
        for smil in smil_list_w_amap:
           scfd_smil = smil
           for idx in range(len(smil_list_wo_amap)):
               lvg_smil = smil_list_wo_amap[idx]
               # check if the scfd and leaving group smiles are the same
               if check_identical_Smiles(scfd_smil,lvg_smil):
                   smil_list_wo_amap[idx] = '<symmetry>'
        # only need the smils that does not have amap
        smil_list_out.append(smil_list_wo_amap)
    return smil_list_out

# get leaving group smiles list and non_amap_total from reaction smiles
#CHANGETO get_RxnSmil2LvgSmil
def get_leaving_group(rxn_smil):
    # separate reactant and product
    rcts,prds=rxn_smil.split('>>')
    rct_list=rcts.split('.')
    # get leaving group from rct using having amap or not in atom
    lvg_list=get_no_AM_group(rct_list)

    non_amap_total=get_noAM_number(prds)
    return lvg_list, non_amap_total


#CHANGETO get_ChangedBond
def show_changed_bond(rct_changed_list_in,prd_changed_list_in):
    changed_bond_list=[]

    # get the changed bond list of rct
    rct_changed_list=[]
    prd_changed_list=[]
    for item in rct_changed_list_in:
        rct_changed_list.extend(item)
    for item in prd_changed_list_in:
        prd_changed_list.extend(item)
    rct_changed_list.sort()
    prd_changed_list.sort()

    # loop through bond_rct and bond_prd to find 
    for bond_rct in rct_changed_list:
        for bond_prd in prd_changed_list:
            not_same_bond_type=True
            # check if the two bonds have same bond type
            for bond_symbol in ['-','=','#',':']:
                if bond_symbol in bond_rct and bond_symbol in bond_prd: 
                    not_same_bond_type = False
            # check if the two bonds with atom map have same bond type as well as the other parts are equal
            if remove_bond_type(bond_rct) == remove_bond_type(bond_prd) and not_same_bond_type:
                # output a bond change expression
                changed_bond_list.append(bond_rct+'>'+bond_prd) 
    # sort the whole list                 
    changed_bond_list.sort()
    return changed_bond_list

# input is list of list of bond
#CHANGETO convert2str_BondListList
def print_str(bond_list_list):
    output= []
    for bond_list in bond_list_list:
        if not bond_list: continue
        output.append(','.join(bond_list))
    return '.'.join(output)

# count number of bond in a bond string expression
#CHANGETO get_BondStr2BondTotal
def count_bond(bond_string):
    count=0
    for character in bond_string:
        if character in  ['-','=','#',':']:
            count+=1
            break
    return count

# get leaving Hs, broken C-H bond, formed C-H bond
#CHANGETO get_Smil2HBondStr
def check_hbond(smil):
    lvg_hs = []
    broken_hbond = []
    formed_hbond = []

    # get rct and prd mol
    allrct, allprd = smil.split('>>')
    rct_list = allrct.split('.')
    prd_list = allprd.split('.')
    allrct_mol = Chem.MolFromSmiles(allrct)
    allprd_mol = Chem.MolFromSmiles(allprd)

    # treat reactant first
    for rct in rct_list:
        lv_h_rct = []
        broken_h = []
        hybrid = get_hybridization(rct)
        rct_mol = Chem.MolFromSmiles(rct)
        # loop through rct and prd to find matching atom by atom map
        for atom1 in rct_mol.GetAtoms():
            for atom2 in allprd_mol.GetAtoms():
                if atom1.GetAtomMapNum() == atom2.GetAtomMapNum():
                    # get hydrogen count
                    atom1_hcount = atom1.GetTotalNumHs()
                    atom2_hcount = atom2.GetTotalNumHs()
                    # if reactant has more hydrogen, means that there is a breaking bond for atom-H
                    if atom1_hcount > atom2_hcount:
                        for i in range(atom1_hcount-atom2_hcount):
                            lv_h_rct.append('H')
                            broken_h.append(atom1.GetSymbol() + hybrid[atom1.GetIdx()]+ '_' + str(atom1.GetAtomMapNum()) + '-H')
        lvg_hs.append(lv_h_rct)
        broken_hbond.append(broken_h)


    # treat reactant first
    for prd in prd_list:
        formed_h = []
        hybrid = get_hybridization(prd)
        prd_mol = Chem.MolFromSmiles(prd)
        # loop through rct and prd to find matching atom by atom map
        for atom1 in prd_mol.GetAtoms():
            for atom2 in allrct_mol.GetAtoms():
                if atom1.GetAtomMapNum() == atom2.GetAtomMapNum():
                    # get hydrogen count
                    atom1_hcount = atom1.GetTotalNumHs()
                    atom2_hcount = atom2.GetTotalNumHs()
                    # if product has more hydrogen, means that there is a forming bond for atom-H
                    if atom1_hcount > atom2_hcount:
                        for i in range(atom1_hcount-atom2_hcount):
                            formed_h.append(atom1.GetSymbol() + hybrid[atom1.GetIdx()]+ '_' + str(atom1.GetAtomMapNum()) + '-H')
        formed_hbond.append(formed_h)
    return lvg_hs, broken_hbond, formed_hbond




#EXAMPLES_FOR_CHECKING____________________________________________________________________________________________

if __name__ == "__main__":
    smil=input('SMILES:')
    #smil='C1[C:2]([OH:1])([c:13]2[cH:14][cH:15][cH:16][cH:17][cH:18]2)[CH2:3]1.O=C1c2ccccc2C(=O)N1O[C:4](=O)[CH2:5][CH2:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1>>[O:1]=[C:2]([CH2:3][CH2:4][CH2:5][CH2:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1)[c:13]1[cH:14][cH:15][cH:16][cH:17][cH:18]1'
    #smil='Br[C:2](=[CH2:1])[CH2:3][N:4]([c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1)[C:11](=[O:12])[C:13]([CH3:14])=[O:15]>>[CH2:1]=[C:2]1[CH2:3][N:4]([c:5]2[cH:6][cH:7][cH:8][cH:9][cH:10]2)[C:11](=[O:12])[C:13]1([CH3:14])[OH:15]'
    #smil='Br[c:17]1[cH:18][cH:19][cH:20][cH:21][c:22]1[N:2]([CH3:1])[C:3](=[O:4])[C:5](=[O:6])[c:7]1[cH:8][cH:9][c:10]2[cH:11][cH:12][cH:13][cH:14][c:15]2[cH:16]1>>[CH3:1][N:2]1[C:3](=[O:4])[C:5]([OH:6])([c:7]2[cH:8][cH:9][c:10]3[cH:11][cH:12][cH:13][cH:14][c:15]3[cH:16]2)[c:17]2[cH:18][cH:19][cH:20][cH:21][c:22]21'
    rct_broken_list, rct_changed_list, prd_formed_list, prd_changed_list = find_diff_bond(smil)
    lv_group, new_group=get_leaving_group(smil)
    changed_bond_list=show_changed_bond(rct_changed_list,prd_changed_list)
    lvg_hs, broken_hbond,formed_hbond = check_hbond(smil)
    for idx_each_rct in range(len(lvg_hs)):
        lv_group[idx_each_rct].extend(lvg_hs[idx_each_rct])
        rct_broken_list[idx_each_rct].extend(broken_hbond[idx_each_rct])
    for idx_each_prd in range(len(formed_hbond)):
        prd_formed_list[idx_each_prd].extend(formed_hbond[idx_each_prd])
    print('broken bonds:', print_str(sort_bond_without_AM(rct_broken_list)))
    print('changing bonds:',sort_bond_without_AM(changed_bond_list,changing_bond=True))
    print('lvg group:', print_str(lv_group))
    print('formed bonds:',print_str(sort_bond_without_AM(prd_formed_list)))
    print('new group:',new_group)