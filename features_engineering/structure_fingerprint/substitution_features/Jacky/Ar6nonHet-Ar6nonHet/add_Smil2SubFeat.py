
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.insert(1,parentdir)
from CORE import *
from get_Info2Feat import *
from get_Smil2SubFeat import *
from get_Smil2Feat import *


#INPUT_FILES_______________________________________________________________________________________________________

integrate_output = True

input_excel_path = '__data_generated/result_16_3_2021_Ar6nonHet-Ar6nonHet_SubINFO_HetINFO.xlsx'
input_target_col_list = [   'OID',
                            'o_Sub_info_rct1',
                            'm_Sub_info_rct1',
                            'p_Sub_info_rct1',
                            'o_Sub_info_rct2',
                            'm_Sub_info_rct2',
                            'p_Sub_info_rct2']

output_excel_path = '__data_generated/result_16_3_2021_Ar6nonHet-Ar6nonHet_SubINFO_HetINFO_SubDescriptor.xlsx'




#FUNCTIONS_______________________________________________________________________________________________________

# break the part of the structure at a distance
def trim_ArSubRingSmil(smil,hassubring_a_id,max_displ = 1):
    output_smil = None

    
    mol = Chem.MolFromSmiles(smil)
    bond_displ_matrix = AllChem.GetDistanceMatrix(mol)
    # the aid connected to * should be = 1
    sub_ac_amap = mol.GetAtomWithIdx(1).GetAtomMapNum()

    # if total_atom_count is less than expected, return original smile
    total_atom_count = len(bond_displ_matrix)-1 # num of atom = all identified atom - 1 from *
    if total_atom_count -1 <= max_displ:    # max_displ count starts from the atom connected to * so -1
        return smil


    # This section is to find the bond to break from substituent atom center (i.e. the target bond break position)
    # get index of neighbor atoms at max_displ
    displ5sub_ac_list = bond_displ_matrix[1]
    # get all atoms at exactly range max_displ (1) from substituent atom center
    start_a_id5sub_ac = [_x for _x in range(len(displ5sub_ac_list)) if displ5sub_ac_list[_x] == max_displ]
    # get all atoms at exactly range max_displ+1 (1+1 = 2) from substituent atom center
    end_a_id5sub_ac = [_x for _x in range(len(displ5sub_ac_list)) if displ5sub_ac_list[_x] == max_displ +1]
    

    # This section is to find the bond to break from the atom that has the same subtituent ring from current atom
    # get index in which radical aid is bond_displ-max_displ away from atom
    bond_displ = bond_displ_matrix[1][hassubring_a_id]
    displ5hassubring_a_list =  bond_displ_matrix[hassubring_a_id]
    # get all atoms at exactly range bond_displ-(max_displ) (becoz the break point locate between the two atoms) from substituent atom center
    start_a_id5hassubring_a = [_x for _x in range(len(displ5hassubring_a_list)) if displ5hassubring_a_list[_x] == bond_displ-max_displ]
    # get all atoms at exactly range bond_displ-(max_displ+1)
    end_a_id5hassubring_a = [_x for _x in range(len(displ5hassubring_a_list)) if displ5hassubring_a_list[_x] == bond_displ-max_displ -1]
    

    # the a_id that satisfy both condition is the only true break point
    start_intxn_a_id = list(set(start_a_id5sub_ac).intersection(start_a_id5hassubring_a))[0]
    end_intxn_a_id = list(set(end_a_id5sub_ac).intersection(end_a_id5hassubring_a))[0]


    # in case if there is an aromatic ring at cutting position, cut the closest nonAromatic bond instead (which should not happen)
    if mol.GetAtomWithIdx(start_intxn_a_id).GetIsAromatic():
        for bond in mol.GetAtomWithIdx(start_intxn_a_id).GetBonds():
            if not bond.GetIsAromatic():
                mol_frag = Chem.FragmentOnBonds(mol,[bond.GetIdx()])
    # normal case: break the bond at the previously identified location
    else:
        mol_frag = Chem.FragmentOnBonds(mol,[mol.GetBondBetweenAtoms(start_intxn_a_id, end_intxn_a_id).GetIdx()])
    
    mol_frag_smil_raw = Chem.MolToSmiles(mol_frag)
    mol_frag_smil_list = mol_frag_smil_raw.split(".")

    # only need the part that has substituent atom center
    for frag_smil in mol_frag_smil_list:
        if ":"+str(sub_ac_amap)+"]" in frag_smil:
            output_smil = frag_smil

    # return error 
    if output_smil == None:
        print("ERROR when trimming ring sub, smiles = ",smil)
        exit()

    return output_smil     



#MAIN_______________________________________________________________________________________________________

if __name__ == '__main__':    
    output_list = []

    # get input from defined excel path
    input_df=pd.read_excel(input_excel_path)
    # extract only needed columns
    input_df_filtered = input_df[input_target_col_list]
    


    electronic_effect_list = []
    steric_effect_list = []

    # all data processing process starts here
    for index, row in input_df_filtered.iterrows():
        rct1_AND_rct2_feat_dict_list = []

        current_oid = row.OID                  #ROWNAME

        rct1_o_sub_smil = row.o_Sub_info_rct1  #ROWNAME
        rct1_m_sub_smil = row.m_Sub_info_rct1  #ROWNAME
        rct1_p_sub_smil = row.p_Sub_info_rct1  #ROWNAME

        rct2_o_sub_smil = row.o_Sub_info_rct2  #ROWNAME
        rct2_m_sub_smil = row.m_Sub_info_rct2  #ROWNAME
        rct2_p_sub_smil = row.p_Sub_info_rct2  #ROWNAME

        rct1_AND_rct2_sub_smil_dict_list =   [  {   "o":rct1_o_sub_smil,
                                                    "m":rct1_m_sub_smil,
                                                    "p":rct1_p_sub_smil,},
                                                {   "o":rct2_o_sub_smil,
                                                    "m":rct2_m_sub_smil,
                                                    "p":rct2_p_sub_smil,}]

        for sub_smil_dict in rct1_AND_rct2_sub_smil_dict_list:
            sub_feat_dict = { "o":[],"m":[], "p":[] }

            # for each position
            for pos in ['o','m','p']:
                sub_smil_raw = sub_smil_dict[pos]
                
                # check if not empty
                if not pd.isna(sub_smil_raw):
                    # split the sub_smil_raw
                    for sub_smil in sub_smil_raw.split("."):
                        
                        # get the info dict of substituent, necessary for identifying RingSub
                        sub_info_dict = get_AllInfo(sub_smil,IsFragment = True,data_type = "smiles")


                        # This part is for treating ring form substituent
                        # get a list of inforamtion about the radical of substiutent
                        # radical is used to label the other ring breaking sites ;assume no radical in original data
                        subring_info_list = [{  'aid'       :_x['aid'],
                                                'radical_e' :_x['radical_e']} for _x in sub_info_dict['atom_info_list']]
                        # get whether there is a ring form substituent
                        hasringsub_boolean = sum([  _x['radical_e'] for _x in subring_info_list])
                        # return a list of radical aid list ()
                        hassubring_a_id_list = [ _x['aid'] for _x in subring_info_list if _x['radical_e'] == 1 ]


                        # get the elec_effect, return None if no info
                        elec_effect = smiles2ArElectronic_Hammett(sub_smil,pos)
                        

                        # get the steric effect
                        # for case has ring substituent
                        if hasringsub_boolean != 0:
                            # assume one radical only
                            hassubring_a_id = hassubring_a_id_list[0]
                            # trim the sub smil to 1 bond displ only, ignore all atoms in the ring after the 1 bond displ (i.e. 2 bond displ from *)
                            sub_smil_trimmed = trim_ArSubRingSmil(sub_smil,hassubring_a_id)

                            # treat the substituent atom center (connecting atom in Ar / *) as a C and account for steric effect
                            sub_smil_trimmed_treated = re.sub("\[\d+\*\]|\(\[\d+\*\]\)","C",sub_smil_trimmed)

                            # error treatment
                            if "*" in sub_smil_trimmed_treated:
                                print("ERROR, there is still * in smiles=",sub_smil_trimmed_treated)
                            
                            # get steric3Dfeat
                            steric3Dfeat = get_Smil2Ster3DFeat(sub_smil_trimmed_treated)

                            # turn off code for StericWholeDescriptor
                            #trimmed_info_dict = {"*":get_AllInfo("[*]"+sub_smil_trimmed,IsFragment = True,data_type = "smiles")}
                            #StericWholeDescriptor = get_Info2SterWholeFeat(trimmed_info_dict)

                        # for case does not have ring substituent
                        else:                            
                            # treat the substituent atom center (connecting atom in Ar / *) as a C and account for steric effect
                            smil_treated = re.sub("\[\d+\*\]|\(\[\d+\*\]\)","C",sub_smil)
                            # if the smiles has * in front
                            smil_treated = re.sub("^\*","C",smil_treated)
                            
                            # get steric3Dfeat
                            steric3Dfeat = get_Smil2Ster3DFeat(smil_treated)

                            # turn off code for StericWholeDescriptor
                            #StericWholeDescriptor = get_Info2SterWholeFeat({"*":sub_info_dict})

                        if elec_effect == None:
                            print("Error: no elec_effect for",sub_smil)
                            exit()

                        sub_feat_dict[pos].append({'elec':elec_effect,'ster':steric3Dfeat})


            # separate the individual feat and how to integrate them tgt to form whole feat for the rct 
            # in case there is some complicated treatment of individual feat in the future
            rct_elec_feat_dict = { "o":0,"m":0, "p":0 }
            rct_ster_feat_dict = { "o":0,"m":0, "p":0 }

            for loc, feat_list in sub_feat_dict.items():
                for feat in feat_list:
                    elec_effect = feat['elec']
                    rct_elec_feat_dict[loc] += elec_effect

                    ster_effect = feat['ster']
                    rct_ster_feat_dict[loc] += ster_effect

            # rct1_o_elec , rct1_m_elec, rct1_p_elec, rct1_o_ster , rct1_m_ster, rct1_p_ster
            rct1_AND_rct2_feat_dict_list+=    [_x for _x in rct_elec_feat_dict.values()] \
                                            + [_x for _x in rct_ster_feat_dict.values()]

        # output = [current_oid, 
        #           rct1_o_electronic,
        #           rct1_m_electronic,
        #           rct1_p_electronic,
        #           rct1_o_steric, 
        #           rct1_m_steric, 
        #           rct1_p_steric,
        #           rct2_o_electronic,
        #           rct2_m_electronic,
        #           rct2_p_electronic,
        #           rct2_o_steric, 
        #           rct2_m_steric, 
        #           rct2_p_steric,]
        
        output = [current_oid] + rct1_AND_rct2_feat_dict_list
        output_list.append(output)


    column_list = [ 'OID',
                    'rct1_o_electronic',
                    'rct1_m_electronic',
                    'rct1_p_electronic',
                    'rct1_o_steric',
                    'rct1_m_steric',
                    'rct1_p_steric',
                    'rct2_o_electronic',
                    'rct2_m_electronic',
                    'rct2_p_electronic',
                    'rct2_o_steric',
                    'rct2_m_steric',
                    'rct2_p_steric']        #ROWNAME

    output_df=pd.DataFrame.from_records(output_list,index=None,columns=column_list)

    if integrate_output == True:
        output_df = pd.merge(input_df, output_df, on=['OID'], how='left')  

    output_df.to_excel(output_excel_path,index=None)

                