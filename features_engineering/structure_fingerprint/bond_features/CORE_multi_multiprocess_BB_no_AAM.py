import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import re
import rdkit.Chem.rdFMCS as MCS
import CORE
from multiprocessing import Process, Manager
from multiprocessing import Pool
from functools import partial

def count_bond(bond_string):
    count=0
    for character in bond_string:
        for bond in ['-','=','#',':']:
            if character==bond:
                count+=1
                break
    return count

def Extract(SMILES):
    #SMILES=input('SMILES:')

    rct_broken_list, rct_changed_list, prd_formed_list, prd_changed_list = CORE.find_diff_bond(SMILES)

    lv_group,NUP,check=CORE.get_leaving_group(SMILES)
    changed_bond_list=CORE.show_changed_bond(rct_changed_list,prd_changed_list)
    leaving_Hs, broken_AH_bond,formed_AH_bond = CORE.check_AH_bond(SMILES)
    if leaving_Hs: #<- check if this list is empty or not
        for idx_each_rct in range(len(leaving_Hs)):
            lv_group[idx_each_rct].extend(leaving_Hs[idx_each_rct])
            rct_broken_list[idx_each_rct].extend(broken_AH_bond[idx_each_rct])
    if formed_AH_bond: #<- check if this list is empty or not
        for idx_each_prd in range(len(formed_AH_bond)):
            try:
                prd_formed_list[0].extend(formed_AH_bond[idx_each_prd]) # 0 because we only do for the case of one product
            except:
                print('ERROR',idx_each_prd,formed_AH_bond,prd_formed_list)
    return rct_broken_list,lv_group,prd_formed_list,changed_bond_list
def cleaning(string,remove_AM=False):
	# Some of this term seem unreasonable, but just a quick fix
    for item in ['Csp1H0','Csp1H1','Csp1H2','Csp1H3','Csp1H4']:
        if item in string:
            string=string.replace(item,'Csp1')
    for item in ['Csp2H0','Csp2H1','Csp2H2','Csp2H3','Csp2H4']:
        if item in string:
            string=string.replace(item,'Csp2')
    for item in ['Csp3H0','Csp3H1','Csp3H2','Csp3H3','Csp3H4']:
        if item in string:
            string=string.replace(item,'Csp3')
    for item in ['Nsp3H0','Nsp3H1','Nsp3H2','Nsp3H3','Nsp3H4']:
        if item in string:
            string=string.replace(item,'Nsp3')
    for item in ['Nsp2H0','Nsp2H1','Nsp2H2','Nsp2H3','Nsp2H4']:
        if item in string:
            string=string.replace(item,'Nsp2')
    #for item in ['F','Br','Cl','I']:
    #    if item in string:
    #        string=string.replace(item,'X')
    if remove_AM:
        string=re.sub(r"_\d+","",string)


    return string


#def do_task(index_list,data_file,breaking_bond,changing_bond,forming_bond,leaving_group,forming_bond_count,NUP_list,check_list,index):
def do_task(output_list,data_file,index):
        row=data_file.iloc[index]
        SMILES=row['AAM']
        rct_broken_list,lv_group,prd_formed_list,changed_bond_list = Extract(SMILES)
        rct_broken_list
        prd_formed_list
        changed_bond_list
        rct_broken_list=cleaning(CORE.print_str(CORE.remove_AM(rct_broken_list)),remove_AM=True)
        prd_formed_list=cleaning(CORE.print_str(CORE.remove_AM(prd_formed_list,remove_CH_bond=False)),remove_AM=True)
        changed_bond_list=cleaning(','.join(CORE.remove_AM(changed_bond_list,changing_bond=True)))


        #print(index, 'done')
        output_list.append([index,rct_broken_list,cleaning(CORE.print_str(lv_group)),changed_bond_list,prd_formed_list,count_bond(prd_formed_list)])
        #output_list.append([index,rct_broken_list])
        print(len(output_list),len(data_file))

def main():
    data_file=pd.read_excel('all_fixed_sofar.xlsx').fillna('')
    output_list = Manager().list()
    input_list=range(0,len(data_file))
    pool = Pool(processes=8)
    # SET INPUT DATA ( NOT INPUT LIST) FOR EACH JOB: RCT_df and PRD_df
    func = partial(do_task,output_list,data_file)
    # SET INPUT LITST FOR EACH, INPUT LIST IS A POOL AND MULTI-PROCESS KEEP PICKING INPUT FROM THIS POOL (edge_list)
    pool.map(func, input_list)
    pool.close()
    pool.join()
    #print(output_list[0])
    extracted_df=pd.DataFrame.from_records(output_list,columns=['index','BREAK_noAAM','LEAVE','CHANGE_noAAM','FORM_noAAM','FORM_count'],index='index')
    extracted_df.sort_index()
    out_df=pd.concat([data_file, extracted_df], axis=1)
    out_df.to_excel('all_fixed_sofar.xlsx',index=None)

if __name__ == '__main__':
    main()
