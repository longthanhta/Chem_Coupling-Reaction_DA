
# Ar 6 member ring (Het or not)
# ring form substituents not considered
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import re

import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.insert(1,parentdir) 

from get_Info import *



#INPUT_FILES_______________________________________________________________________________________________________

label_num = 888

input_excel_path = '__data_raw/result_16_3_2021.xlsx'
input_target_col_list = [   'OID',
                            'SCF_smiles']
                            
output_excel_path = '__data_generated/result_16_3_2021_IsAr6_SubINFO_HetINFO.xlsx'


#FUNCTIONS_______________________________________________________________________________________________________

# Check if reaction center is in defined size of Ar ring 
# the ring_info will combine all ring that has same time as a whole ring system
# therefore, if there is any overlap, the ring must be the opposite type (Ar - nonAr overlap)
def get_ArRingAts(rc_atom,mol,aid2amap_dict,ring_num, IsMapped=True, label_num = label_num):
    # ringat means ring_atoms
    ringats = None
    overlap_w_ring_a_num_list_list = []

    # if reaction center is aromatic
    if rc_atom.GetIsAromatic() == True:
        # get ring info dict for every ring
        ring_info = get_Mol2RingInfo(mol,aid2amap_dict = aid2amap_dict,Mapped=IsMapped)
        for target_ring_info in ring_info:
            # identify the ring that has labelled atom
            if label_num in target_ring_info['ringAts']:
                # check if the ring size is the wanted size
                if target_ring_info['ringSize'] == ring_num:
                    # get all atom id (IsMapped=False) or amap (IsMapped=True)
                    ringats = target_ring_info['ringAts']
                    for other_ring_info in ring_info:
                        overlap_w_ring_a_num_list = []
                        # skip itself
                        if other_ring_info['ringAts'] == ringats:
                            continue
                        
                        # atom_num refers to atom id (IsMapped=False) or amap (IsMapped=True)
                        for atom_num in ringats:
                            # skip target atom
                            if atom_num == label_num:
                                continue
                            
                            # check if it overlaps with other ring (nonAr ring in this case)
                            if atom_num in other_ring_info['ringAts']:
                                overlap_w_ring_a_num_list.append(atom_num)

                        # if there is atom that overlap with other ring (nonAr ring in this case)
                        # ,record the atom_num involved in that ring
                        if overlap_w_ring_a_num_list != []:
                            overlap_w_ring_a_num_list_list.append(overlap_w_ring_a_num_list)
                # does not care about atom other than labelled atom
                break

    # overlap_w_ring_a_num_list_list refers to a list of rings which is expressed by a list of atom num involved in that ring
    return ringats, overlap_w_ring_a_num_list_list

# get normal substituents of any aromatic ring (not consider substituent ring)
def get_ArSub(target_atom,mol,ringats):
    target_a_amap = target_atom.GetAtomMapNum()

    sub_smil_list = []

    for bond in target_atom.GetBonds():
        # get and identify the neighbor atom
        nghbr_atom = bond.GetOtherAtom(target_atom)
        nghbr_a_amap = nghbr_atom.GetAtomMapNum()

        # skip if the neighbor atom is the one in the same ring (Ar ring in this case) i.e. Ar bond
        if nghbr_a_amap in ringats:
            continue

        mol_copy = Chem.Mol(mol)
        # get a sub from cutting bond
        # this method can only cut one bond at a time, thus no ring substituent can be identify
        frag_mol = Chem.FragmentOnBonds(mol_copy,[bond.GetIdx()])
        frag_mol_smil_raw = Chem.MolToSmiles(frag_mol)
        # generate a list of fragment smiles from that cut bond action
        frag_mol_smil_list = frag_mol_smil_raw.split(".")
        
        for frag_mol_smil in frag_mol_smil_list:
            # only skip the fragment with label_num, i.e. only get substituent
            if str(label_num) not in frag_mol_smil:
                sub_smil_list.append(frag_mol_smil)
    
    return (target_a_amap,sub_smil_list)

# only get substituent ring
def get_ArSubRing(target_atom,other_atom_list,mol,ringats):
    target_a_amap = target_atom.GetAtomMapNum()

    sub_smil_list = []
    frag_bond_list = []

    # for all "other_atom" that is in both ring: record the bond to break
    for current_atom in other_atom_list:
        # get and identify the neighbor atom
        current_a_amap = current_atom.GetAtomMapNum()
        for bond in current_atom.GetBonds():
            # get the neighbor atom properties
            nghbr_atom = bond.GetOtherAtom(current_atom)
            nghbr_a_amap = nghbr_atom.GetAtomMapNum()

            # skip if the neighbor atom is the one in the same ring (Ar ring in this case) i.e. Ar bond
            if nghbr_a_amap in ringats:
                continue
                
            # record the bond to break
            frag_bond_list.append(bond.GetIdx())
    
    # cut the bond based on frag_bond_list and unlabel the bond break
    mol_copy = Chem.Mol(mol)
    frag_mol = Chem.FragmentOnBonds(mol_copy,frag_bond_list)
    frag_mol_smil_raw = Chem.MolToSmiles(frag_mol)
    # this is not where we want * to start at so we delete the *
    # remove * and turn the connected C into a radical C
    frag_mol_smil = re.sub("\[\d+\*\]|\(\[\d+\*\]\)","",frag_mol_smil_raw)
    mol2 = Chem.MolFromSmiles(frag_mol_smil)

    # check and stop if any problem occurs
    if mol2 == None:
        print(frag_mol_smil)
        exit()

    # now the previous bond break is unlabelled (does not have a * attached)
    # we cut the bond at target atom to generate a substituent with only one *
    for current_atom in mol2.GetAtoms():
        # get and identify the target atom
        if current_atom.GetAtomMapNum() == target_a_amap:
            for bond in current_atom.GetBonds():
                # get the neighbor atom properties
                nghbr_atom = bond.GetOtherAtom(current_atom)
                nghbr_a_amap = nghbr_atom.GetAtomMapNum()

                # skip if the neighbor atom is the one in the same ring (Ar ring in this case) i.e. Ar bond
                if nghbr_a_amap in ringats:
                    continue

                # cut the bond
                mol3 = Chem.Mol(mol2)
                frag_mol = Chem.FragmentOnBonds(mol3,[bond.GetIdx()])
                frag_mol_smil_raw = Chem.MolToSmiles(frag_mol)
                frag_mol_smil_list = frag_mol_smil_raw.split(".")

    for frag_mol_smil in frag_mol_smil_list:
        # only skip the fragment with label_num, i.e. only get substituent
        if str(label_num) not in frag_mol_smil:
            sub_smil_list.append(frag_mol_smil)
            
    return (target_a_amap,sub_smil_list)



#MAIN_______________________________________________________________________________________________________

if __name__ == '__main__':    
    output_list = []

    # get input from defined excel path
    input_df=pd.read_excel(input_excel_path)
    # extract only needed columns
    input_df_filtered = input_df[input_target_col_list]
    
    # all data processing process starts here
    for index, row in input_df_filtered.iterrows():
        current_oid = row.OID               #ROWNAME
        current_scfd_smil = row.SCF_smiles  #ROWNAME
        rct1_AND_rct2_scfd_smil_list = current_scfd_smil.split(".")
        
        rct1_AND_rct2_IsAr6_boolean_list = []
        rct1_AND_rct2_sub_dict_list = []
        rct1_AND_rct2_het_str_list = []

        # assume there are 2 smils only: one for rct1 another for rct2
        for scfd_smil in rct1_AND_rct2_scfd_smil_list:
            # dedault values for no match
            ismatch_boolean = 0
            rc_atom = None
            het_str = None
            sub_dict = {"o":[],"m":[],"p":[]} 

            # dict for convenience
            #  1 = 1 bond away (i.e. ortho), 2 = 2 bonds away (i.e. meta),3 = 3 bonds away (i.e. para)
            bond_displ2pos_dict = {1:"o", 2:"m", 3:"p"}
            
            # get the mol from smil and start processing
            mol = Chem.MolFromSmiles(scfd_smil)

            # generate a_id a_amap dict for convenience
            aid2amap_dict = {}
            amap2aid_dict = {}
            for atom in mol.GetAtoms():
                a_amap = atom.GetAtomMapNum()
                a_id = atom.GetIdx()
                aid2amap_dict[a_id] = a_amap
                amap2aid_dict[a_amap] = a_id
                if a_amap == label_num:
                    rc_atom = atom
            if rc_atom == None:
                print("ERROR no atom map for",current_scfd_smil)
                exit()

            # get a_amap of all memebrs of Ar ring with size 6 and the substituent ring_info sticking to the Ar6
            # the size cen be changed by ring_num
            ringats, overlap_w_ring_a_num_list_list = get_ArRingAts(rc_atom,mol,aid2amap_dict,ring_num=6,IsMapped=True)
            
            
            if ringats != None:
                # As there is a match of Ar ring
                ismatch_boolean = 1

                # default value
                het_str = []

                # a convenient way to check No. of bond away bewteen any 2 atoms
                bond_displ_matrix = Chem.GetDistanceMatrix(mol)

                for current_a_amap in ringats:
                    # skip if the amap is the reaction center
                    if current_a_amap == label_num:
                        continue

                    # basic info for an atom
                    current_a_id = amap2aid_dict[current_a_amap]
                    current_atom = mol.GetAtomWithIdx(current_a_id)
                    current_a_elem = current_atom.GetSymbol()
                    # find No. of bond away and pos from reaction center
                    bond_displ = int(bond_displ_matrix[amap2aid_dict[label_num]][current_a_id])
                    # convert to current pos from bond displ
                    current_pos = bond_displ2pos_dict[bond_displ]

                    # get Heterocyclic info as a string
                    if current_a_elem != "C":
                        het_str.append(str(bond_displ)+"_"+current_a_elem)

                    # get substituent list of the current atom
                    _ , current_a_sub_list = get_ArSub(current_atom,mol,ringats)
                    # if the normal substituent got is not empty, record it in sub_dict
                    if current_a_sub_list != []:
                        sub_dict[current_pos] += current_a_sub_list

                    # if there is ring substituent
                    for overlap_w_ring_a_num_list in overlap_w_ring_a_num_list_list: 
                        if current_a_amap in overlap_w_ring_a_num_list:
                            other_atom_list = []
                            for other_a_amap in overlap_w_ring_a_num_list:
                                # skip the curre
                                if other_a_amap == current_a_amap:
                                    continue
                                other_a_id = amap2aid_dict[other_a_amap]
                                other_atom = mol.GetAtomWithIdx(other_a_id)
                                other_atom_list.append(other_atom)
                            current_a_sub_ring = get_ArSubRing(current_atom,other_atom_list,mol,ringats)
                            sub_dict[current_pos]+=current_a_sub_ring[1]

                # if no hetero info, return a None
                if het_str == []:
                    het_str = None

            # add data
            rct1_AND_rct2_IsAr6_boolean_list.append(ismatch_boolean)
            rct1_AND_rct2_sub_dict_list.append(sub_dict)
            rct1_AND_rct2_het_str_list.append(het_str)

        # filter applied here
        if rct1_AND_rct2_IsAr6_boolean_list != [0,0]:

            rct1_AND_rct2_sub_info_output_list = []
            for _sub_info in rct1_AND_rct2_sub_dict_list:
                for _, _sub_smil in _sub_info.items():
                    # rct1_o , rct1_m, rct1_p, rct2_o , rct2_m, rct2_p
                    rct1_AND_rct2_sub_info_output_list.append(".".join(_sub_smil))

            # output = [current_oid, 
            #           current_scfd_smil, 
            #           ismatch_boolean (for rct1), 
            #           ismatch_boolean (for rct2),
            #           rct1_o,
            #           rct1_m,
            #           rct1_p,
            #           rct2_o, 
            #           rct2_m, 
            #           rct2_p
            #           het_str (for rct1),
            #           het_str (for rct2),]
            output = [current_oid, current_scfd_smil] + rct1_AND_rct2_IsAr6_boolean_list +rct1_AND_rct2_sub_info_output_list + rct1_AND_rct2_het_str_list
            output_list.append(output)


    column_list = [ 'OID',  
                    'SCF_smiles',
                    'IsAr6_rct1',
                    'IsAr6_rct2',
                    'o_Sub_info_rct1',
                    'm_Sub_info_rct1',
                    'p_Sub_info_rct1',
                    'o_Sub_info_rct2',
                    'm_Sub_info_rct2',
                    'p_Sub_info_rct2',
                    'Het_info_rct1',
                    'Het_info_rct2']        #ROWNAME

    output_df=pd.DataFrame.from_records(output_list,index=None,columns=column_list)

    output_df.to_excel(output_excel_path,index=None)
