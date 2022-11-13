import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolops
import re
# variable list and meaning
# am: atom mapping
# ac: atom Core
# oa: other atom
# b: bond
# nob: number of bond
def position_detection(SCF_smiles_input):

    def ring_double_check(ring_a_lst):
        ring_a_lst_check=[]
        for ra_i in ring_a_lst:
            ra=SCF_mol.GetAtomWithIdx(ra_i)
            cnb_count=0 # correct neightbor count
            #print('atom map',ra.GetAtomMapNum())
            for b in ra.GetBonds():
                oa=rdchem.Bond.GetOtherAtom(b,ra)
                if oa.GetIdx() in ring_a_lst:
                    cnb_count+=1
            if cnb_count>=2:
                #print('correct',cnb_count)
                ring_a_lst_check.append(ra_i)
        ring_a_lst=ring_a_lst_check
        ring_a_lst=list(set(ring_a_lst))
        return ring_a_lst
    # Get Atom Core

    print('SCF_smiles_input_________________________________',SCF_smiles_input)
    SCF_mol=Chem.MolFromSmiles(SCF_smiles_input)

    ac=False
    ac_i=-1 # set -1 because 0 also an index
    for atom in SCF_mol.GetAtoms():
        am=atom.GetAtomMapNum()
        if am == 888: # this value is pre defined in the scaffold extraction script
            ac=atom # ac means atom core
            ac_i=atom.GetIdx()
            break

    # Identify the main ring
    ring_a_lst=[]
    ring_a_lst_1st=[]
    ring_a_lst_2nd=[]
    ring_a_lst_3rd=[]

    ring_a_lst.append(ac_i)
    # Scan three time to get all 6 atoms from reaction center
    # Get index because we cannot trace back atom from atom mapping convincely
    # Get atom map because after cut bond, index may be reset
    for b1 in ac.GetBonds(): # C in ortho positions
        oa1=rdchem.Bond.GetOtherAtom(b1,ac) # oa means other tom
        oa1_archeck=oa1.GetIsAromatic()
        oa1_i=oa1.GetIdx()
        if oa1_archeck and oa1 not in ring_a_lst:
            ring_a_lst_1st.append(oa1_i)
            ring_a_lst.append(oa1_i)
    for ai1 in ring_a_lst_1st: # C in meta positions
        ra1=SCF_mol.GetAtomWithIdx(ai1)
        for b2 in ra1.GetBonds():
            oa2=rdchem.Bond.GetOtherAtom(b2,ra1)
            oa2_archeck=oa2.GetIsAromatic()
            oa2_i=oa2.GetIdx()
            if oa2_archeck and oa2_i not in ring_a_lst:
                ring_a_lst_2nd.append(oa2_i)
                ring_a_lst.append(oa2_i)
    for ai2 in ring_a_lst_2nd: # C in para positions
        ra2=SCF_mol.GetAtomWithIdx(ai2)
        for b3 in ra2.GetBonds():
            oa3=rdchem.Bond.GetOtherAtom(b3,ra2)
            oa3_archeck=oa3.GetIsAromatic()
            oa3_i=oa3.GetIdx()
            if oa3_archeck and oa3_i not in ring_a_lst:
                ring_a_lst_3rd.append(oa3_i)
                ring_a_lst.append(oa3_i)
    print('ring_a_lst_1st',ring_a_lst_1st)
    print('ring_a_lst_1st',ring_a_lst_2nd)
    print('ring_a_lst_1st',ring_a_lst_3rd)

    # Double check and get atom map lst
    # wrong case like: [cH:7]1[cH:8][cH:9][cH:10][c:11](-[c:12]2[n:13][cH:14][cH:15][c:16]3[cH:17][cH:18][cH:19][cH:20][c:21]23)[c:888]1
    # checking by scanning each atom, at least two neightbor in in ring lst
    while len(ring_a_lst)>6:
        ring_a_lst=ring_double_check(ring_a_lst)
    print('ring_a_lst',ring_a_lst)
    ring_a_AM_lst=[]
    # append atom mapping list
    for ra_i in ring_a_lst:
        ring_a_AM_lst.append(SCF_mol.GetAtomWithIdx(ra_i).GetAtomMapNum())
    # Get substitution
    sub_lst=[]
    pos_lst=[]
    for ra_i in ring_a_lst:
        ra=SCF_mol.GetAtomWithIdx(ra_i)
        print('main ring AM',ra.GetAtomMapNum()) # check
        ra_nob=len(ra.GetBonds())
        position = False
        if ra_nob==3:
            # Get the positions
            if ra_i in ring_a_lst_1st: position ='o'
            if ra_i in ring_a_lst_2nd: position ='m'
            if ra_i in ring_a_lst_3rd: position ='p'
            # Get the idex of the atom connect to this caron:
            for b in ra.GetBonds():
                oa=rdchem.Bond.GetOtherAtom(b,ra)
                oa_i=oa.GetIdx()
                oa_m=oa.GetAtomMapNum()
                if oa_i in ring_a_lst: continue # avoid scanning the ring itsel
                # set the atom map of connection atom to 777
                edmol = Chem.EditableMol(SCF_mol)
                edmol.RemoveBond(ra_i, oa_i)
                newmol = edmol.GetMol()
                nm_s_all=Chem.MolToSmiles(newmol) # nm : new mol
                nm_s_lst=nm_s_all.split('.')
            # remove the main ring
                print('nm_s_lst',nm_s_lst)
                for nm_s in nm_s_lst:
                    if '888' in nm_s: continue
                    print('atom index',ra_i,'with position',position)
                    print('has substitution',nm_s)
                    nm_s = re.sub(r'\:\d+', '', nm_s) # remove atom map
                    nm_s = re.sub(r'\:\d+', '', nm_s)
                    print(nm_s)
                    sub_lst.append(nm_s)
                    pos_lst.append(position)
    if ac_i in ring_a_lst_1st: ring_a_lst_1st.remove(ac_i)
    if ac_i in ring_a_lst_2nd: ring_a_lst_2nd.remove(ac_i)
    if ac_i in ring_a_lst_3rd: ring_a_lst_3rd.remove(ac_i)

    sub_s='.'.join(sub_lst)
    pos_s='.'.join(pos_lst)
    return sub_s,pos_s

        #print(no_bond)
if __name__ == '__main__':
    data_df=pd.read_excel('input.xlsx')
    sub_out_lst1=[]
    sub_out_lst2=[]
    pos_out_lst1=[]
    pos_out_lst2=[]
    for index, row in data_df.iterrows():
        #if index!=179:continue # to debug
        print(index,'_______________________________')
        SCF_smiles=row.SCF_smiles
        SCF1,SCF2=SCF_smiles.split('.')
        sub1,pos1=position_detection(SCF1)
        sub2,pos2=position_detection(SCF2)
        sub_out_lst1.append(sub1)
        sub_out_lst2.append(sub2)
        pos_out_lst1.append(pos1)
        pos_out_lst2.append(pos2)
    data_df['sub1']=sub_out_lst1
    data_df['pos1']=pos_out_lst1
    data_df['sub2']=sub_out_lst2
    data_df['pos2']=pos_out_lst2
    data_df.to_excel('out.xlsx',index=None)
