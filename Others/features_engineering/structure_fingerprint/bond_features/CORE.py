import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import re
import rdkit.Chem.rdFMCS as MCS

def Draw_temp(smarts,filename,skip_error=False):
    def Do_Draw_temp(smarts,filename):
        ##print(smarts)
        rxn = Chem.ReactionFromSmarts(smarts)
        rimage = Draw.ReactionToImage(rxn)
        rimage.save(filename)
    if skip_error:
        try:
            Do_Draw_temp(smarts,filename)
        except:
            pass
    else:
        Do_Draw_temp(smarts,filename)

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

def get_hybridization(SMILES):
    element_to_get = ['C','N']
    mol = Chem.MolFromSmiles(SMILES)
    hybrid = ['' for atom in mol.GetAtoms()]
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in element_to_get: continue
        surround_bonds = []
        atom_id = atom.GetIdx()
        for bond in mol.GetBonds():
            aid1 = bond.GetBeginAtomIdx()
            aid2 = bond.GetEndAtomIdx()
            if atom_id == aid1 or atom_id == aid2:
                surround_bonds.append(get_bond_type(bond))
        result = pd.value_counts(surround_bonds)
        if ':' in result.index:
            hybrid[atom_id] = 'sp2'
        elif '#' in result.index:
            hybrid[atom_id] = 'sp1'
        elif '=' in result.index:
            if result['='] == 2:
                hybrid[atom_id] = 'sp1'
            else:
                hybrid[atom_id] = 'sp2'
        else:
            hybrid[atom_id] = 'sp3'
    return hybrid

def get_bond_text(hybrid_atom1,hybrid_atom2,atom_1_map,atom_2_map,element1,element2,no_H1,no_H2,bond_type):
    atom_string_1=element1+hybrid_atom1+'H'+str(no_H1)+'_'+str(atom_1_map)
    atom_string_2=element2+hybrid_atom2+'H'+str(no_H2)+'_'+str(atom_2_map)
    if element1+str(atom_1_map)>element2+str(atom_2_map):
        bond_out=atom_string_1+bond_type+atom_string_2
    else:
        bond_out=atom_string_2+bond_type+atom_string_1
    return bond_out

def get_AM_cut_index(SMILES):
    mol = Chem.MolFromSmiles(SMILES)
    atom_cut=[]
    atom_leave=[]
    for bond in mol.GetBonds():
        aid1 = bond.GetBeginAtomIdx()
        aid2 = bond.GetEndAtomIdx()
        atom1=mol.GetAtomWithIdx(aid1)
        atom2=mol.GetAtomWithIdx(aid2)
        atom_1_map=atom1.GetAtomMapNum()
        atom_2_map=atom2.GetAtomMapNum()
        if atom_1_map==0 and atom_2_map!=0:
            atom_cut.append(aid2)
            atom_leave.append(aid1)
        elif atom_2_map==0 and atom_1_map!=0:
            atom_cut.append(aid1)
            atom_leave.append(aid2)
    return list(set(atom_cut)), list(set(atom_leave))

def get_list_of_bond(SMILES):
    #print('get list of bond SMILES',SMILES)
    bond_list=[] # bond has at least 1 atom map
    unmapped_molecule=True # Check unmapped molecule, if there is a molecule without any mapping atom, ignore it
    mol = Chem.MolFromSmiles(SMILES)
    hybrid = get_hybridization(SMILES)
    for bond in mol.GetBonds():
        #print(bond)
        aid1 = bond.GetBeginAtomIdx()
        aid2 = bond.GetEndAtomIdx()
        bond_type=get_bond_type(bond)
        atom1=mol.GetAtomWithIdx(aid1)
        atom2=mol.GetAtomWithIdx(aid2)
        no_H1=atom1.GetTotalNumHs()
        no_H2=atom2.GetTotalNumHs()
        atom_1_map=atom1.GetAtomMapNum()
        atom_2_map=atom2.GetAtomMapNum()
        element1=atom1.GetSymbol()
        element2=atom2.GetSymbol()
        # if element1 in ['Cl','Br','F','I']:
        #     element1='X'
        # if element2 in ['Cl','Br','F','I']:
        #     element2='X'
        hybrid_atom1 = hybrid[aid1]
        hybrid_atom2 = hybrid[aid2]
        bond_out=get_bond_text(hybrid_atom1,hybrid_atom2,atom_1_map,atom_2_map,element1,element2,no_H1,no_H2,bond_type)
        if not (atom_1_map ==0 and atom_2_map==0):
            unmapped_molecule=False
            bond_list.append(bond_out)
    if unmapped_molecule:
        return []
    return bond_list

def find_diff_bond(reaction):

    rct,prd = reaction.split('>>')
    rct_bond_list=list(get_list_of_bond(rct))
    prd_bond_list=list(get_list_of_bond(prd))


    rct_broken_list = []
    rct_changed_list = []
    prd_formed_list = []
    prd_changed_list = []
    rct_list=[]
    prd_list=[]
    rcts = rct.split('.')
    prds = prd.split('.')

    for rct in rcts:
        bond_list = get_list_of_bond(rct)
        #print('prd_bond_list',prd_bond_list)
        broken_list,changed_list = bond_subtraction(bond_list, prd_bond_list)
        rct_broken_list.append(broken_list)
        rct_changed_list.append(changed_list)
    for prd in prds:
        bond_list = get_list_of_bond(prd)
        formed_list, changed_list = bond_subtraction(bond_list, rct_bond_list)
        prd_formed_list.append(formed_list)
        prd_changed_list.append(changed_list)


    return rct_broken_list, rct_changed_list, prd_formed_list, prd_changed_list

def remove_AM_list(bond_list):
    bond_out=[]
    for bond in bond_list:
        bond_out.append(re.sub(r"_\d+","",bond))
    return bond_out
def remove_AM_single(bond):
    bond_no_AM=[]
    return re.sub(r"_\d+","",bond)

def sort_bond(bond):
    for connect in '-=#:':
        if connect in bond:
            break
    c = bond.index(connect)
    a1 = bond[:c]
    a2 = bond[c+1:]
    [a1,a2] = sorted([a1,a2])
    return a1+connect+a2


def remove_AM(input_list,list_of_list=True,changing_bond=False,remove_CH_bond=False,with_AM=False):
    if not list_of_list:
        bond_list_out=[]
        for bond in input_list:
            if remove_CH_bond:
                if '-H' in bond or 'H-' in bond:
                    continue
            if not with_AM:
               #bond_list_out.append(sort_bond(re.sub(r"_\d+","",bond)))
               bond_list_out.append(bond)

            else:
               bond_list_out.append(bond)
        bond_list_out
        return bond_list_out
    if changing_bond:
        bond_list_out=[]
        for couple_bond in input_list:
            bond_list=couple_bond.split('>')
            bond_list_out.append(re.sub(r"_\d+","",bond_list[0])+'>'+re.sub(r"_\d+","",bond_list[1]))
        bond_list_out
        return bond_list_out
    else:
        bond_list_out=[]
        for mol in input_list:
            mol_out=[]
            for bond in mol:
                if remove_CH_bond:
                    if '-H' in bond or 'HH0-' in bond:
                        continue
                #mol_out.append(sort_bond(re.sub(r"_\d+","",bond)))
                mol_out.append(bond)
            mol_out
            bond_list_out.append(mol_out)
        bond_list_out
        return bond_list_out
def remove_bond_type(bond):
    bond_no_bt=[]
    for bond_symbol in ['-','=','#',':']:
        if bond_symbol in bond:
            bond_type = bond_symbol
    atoms = bond.split(bond_type)
    atoms = [atom.split('_')[0].split('H')[0].split('sp')[0]+'_'+atom.split('_')[1] for atom in atoms]
    atoms
    bond_no_bt='-'.join(atoms)
    return bond_no_bt

def bond_subtraction(list1, list2,debug_mode=False): # Set False for debug mode
    different_list = list(set(list1) - set(list2))
    remove_bt_list2 = [remove_bond_type(bond) for bond in list2]
    new = []
    changed = []
    for bond in different_list:
        if remove_bond_type(bond) in remove_bt_list2:
            changed.append(bond)
        else:
            new.append(bond)
    if debug_mode:
            new = remove_AM_list(new)
    return new, changed


def get_AM_number(SMILES):
    cnt = 0
    mol = Chem.MolFromSmiles(SMILES)
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() != 0:
            cnt += 1
    return cnt

def get_noAM_number(SMILES):
    cnt = 0
    mol = Chem.MolFromSmiles(SMILES)
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            cnt += 1
    return cnt

def check_identical_Smiles(SMILES1, SMILES2):
    mol1 = Chem.MolFromSmiles(SMILES1)
    mol2 = Chem.MolFromSmiles(SMILES2)
    try:
        return mol2.HasSubstructMatch(mol1) and mol1.HasSubstructMatch(mol2)
    except:
        return False

def get_no_AM_group(SMILES_list):
    SMILES_list_out = []
    for each_SMILES in SMILES_list:
        #print('RCT',RCT)
        atom_cut, atom_leave = get_AM_cut_index(each_SMILES)
        SMILES_list_no_AM=[]
        SMILES_list_with_AM=[]
        #print('RCT',RCT)
        for leave_point in atom_leave:
            mw=Chem.RWMol(Chem.MolFromSmiles(each_SMILES))
            mw.RemoveAtom(leave_point)
            SMILES_list=Chem.MolToSmiles(mw).split('.')
            for SMILES in SMILES_list:
            #print('SMILES',SMILES)
                contain_AM=False
                if re.search(':\d+', SMILES):
                    contain_AM=True
                if contain_AM:
                    SMILES_list_with_AM.append(SMILES)
        for cut_point in atom_cut:
            mw=Chem.RWMol(Chem.MolFromSmiles(each_SMILES))
            mw.RemoveAtom(cut_point)
            SMILES_list=Chem.MolToSmiles(mw).split('.')
            for SMILES in SMILES_list:
                #print('SMILES',SMILES)
                contain_AM=False
                if re.search(':\d+', SMILES):
                    contain_AM=True
                if not contain_AM:
                    SMILES_list_no_AM.append(SMILES)
        #print(SMILES_list_with_AM)
        for SMILES in SMILES_list_with_AM:
           staying = SMILES
           for idx in range(len(SMILES_list_no_AM)):
               leaving = SMILES_list_no_AM[idx]
               #if check_identical_Smiles(staying,leaving):
               #    SMILES_list_no_AM[idx] = '[symmetry]'
        SMILES_list_out.append(SMILES_list_no_AM)

    return SMILES_list_out

def get_leaving_group(reaction_SMILES):
    RCTs,PRDs=reaction_SMILES.split('>>')
    RCT_list=RCTs.split('.')
    RCT_list_out=get_no_AM_group(RCT_list)
    no_AM_PRD=get_noAM_number(PRDs)
    AM_PRD=get_AM_number(PRDs)
    AM_RCT=get_AM_number(RCTs)
    check = (AM_PRD>AM_RCT)
    return RCT_list_out, no_AM_PRD, check


def show_changed_bond(rct_changed_list_in,prd_changed_list_in):
    changed_bond_list=[]
    rct_changed_list=[]
    prd_changed_list=[]
    for item in rct_changed_list_in:
        rct_changed_list.extend(item)
    for item in prd_changed_list_in:
        prd_changed_list.extend(item)
    rct_changed_list
    prd_changed_list
    for bond_r in rct_changed_list:
        for bond_p in prd_changed_list:
            not_same_bond=True
            for bond_symbol in ['-','=','#',':']:
                if bond_symbol in bond_r and bond_symbol in bond_p:
                    not_same_bond = False
            if remove_bond_type(bond_r) == remove_bond_type(bond_p) and not_same_bond:
                changed_bond_list.append(bond_r+'>'+bond_p)
    changed_bond_list
    return changed_bond_list

def print_str(bond_list):
    output= []
    for bonds in bond_list:
        if not bonds: continue
        output.append(','.join(bonds))
    return '.'.join(output)
def count_bond(bond_string):
    count=0
    for character in bond_string:
        for bond in ['-','=','#',':']:
            if character==bond:
                count+=1
                break
    return count
def check_AH_bond(SMILES):
    leaving_Hs = []
    broken_AH_bond = []
    formed_AH_bond = []
    rct, prd = SMILES.split('>>')
    rcts = rct.split('.')
    prds = prd.split('.')
    rctmol = Chem.MolFromSmiles(rct)
    prdmol = Chem.MolFromSmiles(prd)
    for rct in rcts:
        lv_H_each_rct = []
        broken_AH_each_rct = []
        formed_AH_each_rct=[]
        hybrid = get_hybridization(rct)
        rctmol = Chem.MolFromSmiles(rct)
        for atom1 in rctmol.GetAtoms():
            for atom2 in prdmol.GetAtoms():
                if atom1.GetAtomMapNum() == atom2.GetAtomMapNum():
                    H1 = atom1.GetTotalNumHs()
                    H2 = atom2.GetTotalNumHs()
                    if H1 > H2:
                        for i in range(H1-H2):
                            lv_H_each_rct.extend('H')
                            broken_AH_each_rct.append(atom1.GetSymbol() + hybrid[atom1.GetIdx()]+ '_' + str(atom1.GetAtomMapNum()) + '-H')
                    if H2 > H1:
                        for i in range(H2-H1):
                            formed_AH_each_rct.append(atom1.GetSymbol() + hybrid[atom1.GetIdx()]+ '_' + str(atom1.GetAtomMapNum()) + '-H')
        leaving_Hs.append(lv_H_each_rct)
        broken_AH_bond.append(broken_AH_each_rct)
        formed_AH_bond.append(formed_AH_each_rct)
        formed_AH_bond=[b for b in formed_AH_bond if b] #<= rmeove empty formed AH-bond, result in a loop of 2 reactants
    if len(broken_AH_bond)==0: broken_AH_bond=False
    if len(formed_AH_bond)==0: formed_AH_bond=False
    return leaving_Hs, broken_AH_bond, formed_AH_bond


if __name__ == "__main__":
    #SMILES=input('SMILES:')
    SMILES='[CH2:1]=[CH:2][c:3]1[cH:12][cH:11][c:6]([C:7]([O:9][CH3:10])=[O:8])[cH:5][cH:4]1.[CH3:13][O:14][C:15]([c:17]([c:18]([H])[n:19]2[CH3:20])[c:26]3[c:21]2[cH:22][cH:23][cH:24][cH:25]3)=[O:16]>>[CH3:10][O:9][C:7]([c:6]4[cH:11][cH:12][c:3]([CH:2]([c:18]([c:17]5[C:15]([O:14][CH3:13])=[O:16])[n:19]([CH3:20])[c:21]6[c:26]5[cH:25][cH:24][cH:23][cH:22]6)[CH2:1][H])[cH:4][cH:5]4)=[O:8]'
    rct_broken_list, rct_changed_list, prd_formed_list, prd_changed_list = find_diff_bond(SMILES)
    lv_group, new_group, check=get_leaving_group(SMILES)
    changed_bond_list=show_changed_bond(rct_changed_list,prd_changed_list)
    leaving_Hs, broken_AH_bond,formed_AH_bond = check_AH_bond(SMILES)
    print('formed_AH_bond',formed_AH_bond)
    print('leaving_Hs',leaving_Hs)
    print('broken_AH_bond',broken_AH_bond)

    for idx_each_rct in range(len(leaving_Hs)):
        lv_group[idx_each_rct].extend(leaving_Hs[idx_each_rct])
        rct_broken_list[idx_each_rct].extend(broken_AH_bond[idx_each_rct])
    for idx_each_prd in range(len(formed_AH_bond)):
        prd_formed_list[idx_each_prd].extend(formed_AH_bond[idx_each_prd])
    print('broken bonds:', print_str(remove_AM(rct_broken_list)))
    print('changing bonds:',remove_AM(changed_bond_list,changing_bond=True))
    print('leaving group:', print_str(lv_group))
    print('formed bonds:',print_str(prd_formed_list))
    print('new group:',new_group)
    print('remember to remove H in the excel file of forming bond')
