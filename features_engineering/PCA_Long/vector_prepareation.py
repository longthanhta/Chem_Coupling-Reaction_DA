import pandas as pd
import pickle
import numpy as np
import rdkit
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import ChiralType
from IPython.display import SVG
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw

def Get_MFP(smiles,nBits=2048):
    fp=(np.array(AllChem.GetMorganFingerprintAsBitVect(
    AllChem.MolFromSmiles(smiles), 1
    , nBits,useChirality=False), dtype=np.int))
    fp=np.array(fp)
    lst=[]
    for i in range(0,len(fp)):
        if fp[i] == 1:
            lst.append(i+1)
    lst.append(smiles)
    return fp


csp1_list=['Alkyne','C#N']
csp2_list=['Vinyl','C=O','Vinylene',]
csp2_aromatic_list=['Het','Ar']

def countringar(SCFstring):
    total_ring=0
    for i in SCFstring:
        if i == 'Ar':
            total_ring = total_ring + 1
    return total_ring

def countC(smiles):
    total_C=0
    for i in smiles:
        if i == 'C' or i == 'c':
            total_C = total_C + 1
    return total_C

def countringhet(SCFstring):
    total_ring=0
    for i in SCFstring:
        if i == 'Het':
            total_ring = total_ring + 1
    return total_ring


data=pd.read_excel('result_classification.xlsx')
print(data)
EGC_lst=[]
SCF_lst=[]
SCFA_lst=[]
# Get unique data
for index,row in data.iterrows():
    EGC1,EGC2=row['EG_class'].split('.')
    EGC_lst.append(EGC1)
    EGC_lst.append(EGC2)
    SCF1,SCF2=row['SCF_gen_uni'].split('.')
    SCF_lst.append(SCF1)
    SCF_lst.append(SCF2)
    #print(row['FORM_noAAM'])
    SCFA1,SCFA2=row['FORM_noAAM'].split('-')
    SCFA_lst.append(SCFA1)
    SCFA_lst.append(SCFA2)
EGC_lst=list(set(EGC_lst))
SCF_lst=list(set(SCF_lst))
SCFA_lst=list(set(SCFA_lst))
SCF_hbdz_lst=['Csp1','Csp2','Csp2Ar','Csp3']
LIG_lst=['P;P;P','P-P','P-P-P']

print(EGC_lst)
print(SCF_lst)
print(SCFA_lst)
EGC_lst.sort()
SCF_lst.sort()
SCFA_lst.sort()
#pd.DataFrame(EGC_lst).to_csv('EGC_lst.csv',index=None,header=['Item'])
#pd.DataFrame(SCF_lst).to_csv('SCF_lst.csv',index=None,header=['Item'])
#pd.DataFrame(SCFA_lst).to_csv('SCFA_lst.csv',index=None,header=['Item'])
EGC_En=pd.read_csv('EGC_lst.csv')
SCF_ring=pd.read_csv('SCF_lst.csv')

# Get bit vector
bit_total_lst=[]
bit_total_MGFP_lst=[]
for index,row in data.iterrows():
    EGC_bit = [0] * len(EGC_lst)
    SCF_bit = [0] * len(SCF_lst)
    SCFA_bit = [0] * len(SCFA_lst)
    SCF_hbdz_bit=[0] * len(SCF_hbdz_lst)
    LIG_bit=[0]*len(LIG_lst)
    RCT_smiles=row['RCT']
    EGC1,EGC2=row['EG_class'].split('.')
    SCF1,SCF2=row['SCF_gen_uni'].split('.')
    SCFA1,SCFA2=row['FORM_noAAM'].split('-')
    LIG=row['LIG_GROUP']
    LIG_smiles=row['LIG_SMILES']
    LIG_type='mono'
    if LIG=='bi': LIG_type='bi'
    if LIG=='tri': LIG_type='tri'

    #print(EGC1,EGC2,SCF1,SCF2,SCFA1,SCFA2)
    EGC_bit[EGC_lst.index(EGC1)]=1
    SCF_bit[SCF_lst.index(SCF1)]=1
    SCFA_bit[SCFA_lst.index(SCFA1)]=1
    EGC_bit[EGC_lst.index(EGC2)]=1
    SCF_bit[SCF_lst.index(SCF2)]=1
    SCFA_bit[SCFA_lst.index(SCFA2)]=1
    LIG_bit[LIG_lst.index(LIG)]=1
    EGC_vec=[]
    for id,it in enumerate(EGC_bit):
        selected_value=EGC_En['en'][EGC_En['Item']==EGC_lst[id]]
        #print('debug',EGC_lst[id],float(selected_value))
        bit=it*selected_value
        bit=bit.tolist()[0]
        #print('bit timed en',bit)
        EGC_vec.append(bit)
    print('EGC_vec',EGC_vec)
    #for id,it in enumerate(SCF_bit):
    #    SCF_vec.append(it*SCF_ring['ring'][id])


    hpbdz1='Csp3'
    if SCF1 in csp1_list:
        hpbdz1='Csp1'
    if SCF1 in csp2_list:
        hpbdz1='Csp2'
    if SCF1 in csp2_aromatic_list:
        hpbdz1='Csp2Ar'
    hpbdz2='Csp3'
    if SCF2 in csp1_list:
        hpbdz2='Csp1'
    if SCF2 in csp2_list:
        hpbdz2='Csp2'
    if SCF2 in csp2_aromatic_list:
        hpbdz2='Csp2Ar'
    print('hpbdz1',hpbdz1)
    print('hpbdz2',hpbdz2)
    SCF_hbdz_bit[SCF_hbdz_lst.index(hpbdz1)]=1
    SCF_hbdz_bit[SCF_hbdz_lst.index(hpbdz2)]=1


    SCF_hbdz_vec=[]
    print('SCF1',SCF1)
    SCF_string=SCF1+'.'+SCF2
    number_of_ring=SCF_string.count('Ar')+SCF_string.count('Het')
    print('number of ring',number_of_ring)
    for id,it in enumerate(SCF_hbdz_bit):
        SCF_hbdz_vec.append(it*number_of_ring)
    print('LIG_smiles',LIG_smiles)
    ligand_weight=LIG_smiles.count('c')+LIG_smiles.count('C')
    LIG_vec=[]
    for id,it in enumerate(LIG_bit):
        LIG_vec.append(it*ligand_weight)
    print(len(SCF_hbdz_bit),'bit hbdz',SCF_hbdz_bit)
    #MGFB_bit=Get_MFP(RCT_smiles)
    #bit_total_MGFP=EGC_bit+SCF_bit+SCFA_bit+MGFB_bit[0]
    #bit_total_MGFP_lst.append(bit_total_MGFP)
    #bit_total=EGC_bit+SCF_bit+SCF_hbdz_bit
    #bit_total=EGC_vec+SCF_bit+SCF_hbdz_vec+LIG_vec
    bit_total=EGC_vec+SCF_hbdz_vec+LIG_vec
    bit_total_lst.append(bit_total)
    #bit_total_MGFP_lst.append(bit_total_MGFP)
    print(len(bit_total),'bit vector',bit_total)

    #break
#header_lst=EGC_lst+SCF_lst+SCF_hbdz_lst+LIG_lst
header_lst=EGC_lst+SCF_hbdz_lst+LIG_lst

# LG:
anotation_lst=[]
for index,row in data.iterrows():
    EG1,EG2=row['EG_class'].split('.')
    EG_lst=[EG1,EG2]
    EG_lst.sort()
    EG_string='.'.join(EG_lst)
    anotation_lst.append(EG_string)
print('headers',len(header_lst),header_lst)
df_out=pd.DataFrame.from_records(bit_total_lst,index=None,columns=header_lst)
df_out['class']=data['classfication']
df_out['EGC']=anotation_lst
df_out['SCF']=data['SCF_gen_uni']
df_out['SCFA']=data['FORM_noAAM']
df_out['LIG']=data['LIG']
df_out['LIG_GROUP']=data['LIG_GROUP']
df_out['RXN']=data['RXN']
df_out['yield']=data['Yield']
df_out['ID']=data['OID']
df_out['TEXT']=anotation_lst
pd.DataFrame(header_lst).to_csv('header_lst.csv',index=None,header=['Item'])

with open('bits.pkl', 'wb') as f:
    pickle.dump(bit_total_lst, f)
#with open('bitsMG.pkl', 'wb') as f:
#    pickle.dump(bit_total_MGFP_lst, f)
print('df_out',df_out)
df_out.to_excel('out.xlsx',index=None)
