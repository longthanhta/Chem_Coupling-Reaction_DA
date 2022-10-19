from rdkit.Chem import AllChem as Chem
import pandas as pd

# putting smarts here
smarts_input='[c:8]:[c:7]:[c;H0;D3;+0:6](:[c:9]:[c:10])-[c;H0;D3;+0:1](:[c:2]:[c:3]):[c:4]:[c:5]>>F-C(=O)-[c;H0;D3;+0:1](:[c:2]:[c:3]):[c:4]:[c:5].O-B(-O)-[c;H0;D3;+0:6](:[c:7]:[c:8]):[c:9]:[c:10]'
rxn = Chem.ReactionFromSmarts(smarts_input)

data_df=pd.read_excel('prd_lst.xlsx')
rxn_smi_lst1=[]
rxn_smi_lst2=[]
prd_lst=[]
for index,row in data_df.iterrows():
    prd_smi=row['prd']
    #print('prd_smi',prd_smi)
    target=prd_smi
    ps = rxn.RunReactants((Chem.MolFromSmiles(target),))
    ps=list(set(ps))
    #print('ps',ps)
    rxn_smi_1case=[]
    for i in range(0,len(ps)):
        #print('case',i)
        rct1_smi=Chem.MolToSmiles(ps[i][0])
        rct2_smi=Chem.MolToSmiles(ps[i][1])
        rxn_smi=rct1_smi+'.'+rct2_smi+'>>'+prd_smi
        rxn_smi_1case.append(rxn_smi)
    rxn_smi_1case=list(set(rxn_smi_1case))
    print('rxn_smi_1case',rxn_smi_1case)
    if len(rxn_smi_1case)==2:
        rxn_smi_lst1.append(rxn_smi_1case[0])
        rxn_smi_lst2.append(rxn_smi_1case[1])
    if len(rxn_smi_1case)==1:
        rxn_smi_lst1.append(rxn_smi_1case[0])
        rxn_smi_lst2.append('na')
data_df['rxn1']=rxn_smi_lst1
data_df['rxn2']=rxn_smi_lst2
data_df.to_excel('rxn_lst.xlsx',index=None)
print('output file is rxn_lst.xlsx')
print('For aryl aryl case, there can be two case since the smart treat both subtrate the same')
