# This script is for extractin scaffold
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import re
import rdkit.Chem.rdFMCS as MCS
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit import DataStructs
import numpy as np
from rdkit.Chem import rdchem
import regex as re
from random import randint

def get_en(elmt_check):
    if elmt_check=='H': return	2.2
    if elmt_check=='Li': return	0.98
    if elmt_check=='Be': return	1.57
    if elmt_check=='B': return	2.04
    if elmt_check=='C': return	2.55
    if elmt_check=='N': return	3.04
    if elmt_check=='O': return	3.44
    if elmt_check=='F': return	3.98
    if elmt_check=='Na': return	0.93
    if elmt_check=='Mg': return	1.31
    if elmt_check=='Al': return	1.61
    if elmt_check=='Si': return	1.9
    if elmt_check=='P': return	2.19
    if elmt_check=='S': return	2.58
    if elmt_check=='Cl': return	3.16
    if elmt_check=='K': return	0.82
    if elmt_check=='Ca': return	1
    if elmt_check=='Sc': return	1.36
    if elmt_check=='Ti': return	1.54
    if elmt_check=='V': return	1.63
    if elmt_check=='Cr': return	1.66
    if elmt_check=='Mn': return	1.55
    if elmt_check=='Fe': return	1.83
    if elmt_check=='Co': return	1.88
    if elmt_check=='Ni': return	1.91
    if elmt_check=='Cu': return	1.9
    if elmt_check=='Zn': return	1.65
    if elmt_check=='Ga': return	1.81
    if elmt_check=='Ge': return	2.01
    if elmt_check=='As': return	2.18
    if elmt_check=='Se': return	2.55
    if elmt_check=='Br': return	2.96
    if elmt_check=='Kr': return	3
    if elmt_check=='Rb': return	0.82
    if elmt_check=='Sr': return	0.95
    if elmt_check=='Y': return	1.22
    if elmt_check=='Zr': return	1.33
    if elmt_check=='Nb': return	1.6
    if elmt_check=='Mo': return	2.16
    if elmt_check=='Tc': return	1.9
    if elmt_check=='Ru': return	2.2
    if elmt_check=='Rh': return	2.28
    if elmt_check=='Pd': return	2.2
    if elmt_check=='Ag': return	1.93
    if elmt_check=='Cd': return	1.69
    if elmt_check=='In': return	1.78
    if elmt_check=='Sn': return	1.96
    if elmt_check=='Sb': return	2.05
    if elmt_check=='Te': return	2.1
    if elmt_check=='I': return	2.66
    if elmt_check=='Xe': return	2.6
    if elmt_check=='Cs': return	0.79
    if elmt_check=='Ba': return	0.89
    if elmt_check=='La': return	1.1
    if elmt_check=='Ce': return	1.12
    if elmt_check=='Pr': return	1.13
    if elmt_check=='Nd': return	1.14
    if elmt_check=='Sm': return	1.17
    if elmt_check=='Gd': return	1.2
    if elmt_check=='Dy': return	1.22
    if elmt_check=='Ho': return	1.23
    if elmt_check=='Er': return	1.24
    if elmt_check=='Tm': return	1.25
    if elmt_check=='Lu': return	1.27
    if elmt_check=='Hf': return	1.3
    if elmt_check=='Ta': return	1.5
    if elmt_check=='W': return	2.36
    if elmt_check=='Re': return	1.9
    if elmt_check=='Os': return	2.2
    if elmt_check=='Ir': return	2.2
    if elmt_check=='Pt': return	2.28
    if elmt_check=='Au': return	2.54
    if elmt_check=='Hg': return	2
    if elmt_check=='Tl': return	1.62
    if elmt_check=='Pb': return	2.33
    if elmt_check=='Bi': return	2.02
    if elmt_check=='Po': return	2
    if elmt_check=='At': return	2.2
    if elmt_check=='Ra': return	0.9
    if elmt_check=='Ac': return	1.1
    if elmt_check=='Th': return	1.3
    if elmt_check=='Pa': return	1.5
    if elmt_check=='U': return	1.38
    if elmt_check=='Np': return	1.36
    if elmt_check=='Pu': return	1.28
    if elmt_check=='Am': return	1.3
    if elmt_check=='Cm': return	1.3
    if elmt_check=='Bk': return	1.3
    if elmt_check=='Cf': return	1.3
    if elmt_check=='Es': return	1.3
    if elmt_check=='Fm': return	1.3
    if elmt_check=='Md': return	1.3
    if elmt_check=='No': return	1.3

#Get a system of rings that have the common bond
def GetSystem(smiles, atom):
    mol = Chem.MolFromSmiles(smiles)
    ri = mol.GetRingInfo()
    systems = []
    idx = atom.GetIdx()
    for ring in ri.AtomRings():
        if not isRingAromatic(mol, ring):
            continue
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon:
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    #print(systems)
    for system in systems:
        if idx in system:
            return system
    return {}


def isRingAromatic(mol, atomRing):
    for id in atomRing:
        if not mol.GetAtomWithIdx(id).GetIsAromatic():
            return False
    return True

#count num of single rings in a system and num of hetero atom then return type
def detect(smiles, atom):
    mol = Chem.MolFromSmiles(smiles)
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    system = GetSystem(smiles, atom)
    #print(system)
    numHet = 0
    numring = 0
    Het = 0
    for ring in rings:
        if not isRingAromatic(mol, ring):
            continue
        ringAts = set(ring)
        if len(ringAts.intersection(system)):
            #print("yes")
            numring +=1
            numHet = 0
            for i in ring:
                if mol.GetAtomWithIdx(i).GetSymbol() != 'C':
                    numHet+=1
                    break
            if numHet > 0:
                Het+=1

    if Het > 0:
        #if Het == 1 and numring == 1:
            #return 'Het'
        #if numHet > 1 and numring > 1:
            #return 'Het2'
        #else:
            #return 'ArHet'
        #Het is Het only
        return 'Het'
    else:
        if numring == 3:
            return 'Ar3'
        if numring == 2:
            return 'Ar2'
        else:
            return 'Ar'


#check type, if C-Aryl, then check whether hetero or multiring
def GetType(smiles, atom, otherAtom):
    # Detect the type
    print('GetType input atom:',atom.GetSymbol())
    print('GetType input other atom',otherAtom.GetSymbol())
    #print('hybridization:',atom.GetHybridization())
    if atom.GetIsAromatic(): # -> C is aromatic -> return
        #print('DEBBUG:Gettype: atom is aromatic')
        return detect(smiles, atom)
    # from now is if C not aromatic
    bonds = atom.GetBonds()
    numAr = 0
    for bond in bonds: # this is bond of the atom
        #print('bond',bond.GetBondType())
        #print(bond.GetBondType() == rdchem.BondType.DOUBLE)
        otherAtom = bond.GetOtherAtom(atom)
        sym_oat = otherAtom.GetSymbol()
        noh_oat = otherAtom.GetTotalNumHs()
        aam_oat=otherAtom.GetAtomMapNum()
        if aam_oat==0: continue
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            if sym_oat == 'C':
                #print('DEBBUG:Gettype: Vinyl')
                bonds_oat=otherAtom.GetBonds()
                if len(bonds_oat)==2:
                    return 'Vinyl'
                else:
                    return 'Vinylene'
            else:
                #print('DEBBUG:Gettype: C=+...')
                return 'C=' + sym_oat
        elif bond.GetBondType() == rdchem.BondType.TRIPLE:
            if sym_oat == 'C':
                #print('DEBBUG:Gettype: Alkyne')
                return 'Alkyne'
            else:
                #print('DEBBUG:Gettype: Alkyne+...')
                return 'C#' + sym_oat
        elif bond.GetBondType() == rdchem.BondType.SINGLE: # Get the alpha arylation case
            if sym_oat =='C':
                bonds_oat=otherAtom.GetBonds()
                for bondoa in bonds_oat:
                    # other of other atom
                    OOA=bondoa.GetOtherAtom(otherAtom)
                    if OOA.GetSymbol()=='O' and bondoa.GetBondType() == rdchem.BondType.DOUBLE:
                        return 'C-(C=O)'
        if otherAtom.GetIsAromatic():
            numAr+=1
    if numAr > 0:
        #print('DEBBUG:Gettype: C+Aryl')
        return 'C-' + 'Aryl'
    #print('DEBBUG: Detect nospecial case-> alkyl')
    return 'Alkyl'

def checkAllyl(at, idx):
    neighbors = at.GetNeighbors()
    for atom in neighbors:
        if atom.GetIdx() == idx:
            continue
        if atom.GetSymbol() == 'C':
            bonds = atom.GetBonds()
            for bond in bonds:
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    if rdchem.Bond.GetOtherAtom(bond, atom).GetSymbol() == 'C':
                        return True
    return False
def string_swap(string):
    CP1,CP2=string.split('.')
    return CP2+'.'+CP1
def get_main(mol, idx1, idx2, otherAtom):
    # Run through each cutting atom
    at = mol.GetAtomWithIdx(idx1)
    oat=otherAtom
    typ = GetType(smiles, at, oat)
    for atom in mol.GetAtoms():
        if atom.GetIdx()==idx1:
            atom.SetAtomMapNum='888'
        else:
            atom.SetAtomMapNum=0
    if typ == 'Alkyl' and checkAllyl(at, idx2):
        typ = 'Allyl'
    if typ == 'C-Aryl':
        t = []
        neighbors = at.GetNeighbors()
        for atom in neighbors:
            if atom.GetIdx() == idx2:
                continue
            if atom.GetIsAromatic():
                t.append(detect(smiles, atom))
        typ = 'C-'
        for i in t:
            typ = typ + '(' + i + ')'
    return typ,Chem.MolToSmiles(mol)

def Get_elem(s):
    if '_' not in s:
        return 'H'
    if ('C' in s or 'Aryl' in s) and 'Cl' not in s:
        return 'C'
    if 'N' in s:
        return 'N'
    return s.split('_')[0]

def Get_emap(s):
    if '_' not in s:
        return 0
    return int(s.split('_')[1])


def get_lg_n_scf_smiles(mol, idx1, idx2):
    print('Start debugging get_lg_n_scf_smiles')
    print('index 2 belogn to other atom:',idx2)
    for atom in mol.GetAtoms():
        if atom.GetIdx()==idx1:
            print('atom on scaffold side atom map no:',atom.GetAtomMapNum())
            atom.SetAtomMapNum(888) # assign index 999 to the atom oflG
        elif atom.GetIdx()==idx2:
            print('atom on leavingg group side atom map no:',atom.GetAtomMapNum())
            atom.SetAtomMapNum(999) # assign index 999 to the atom oflG
    edmol = Chem.EditableMol(mol)
    if idx2>-1:
        edmol.RemoveBond(idx1, idx2)
    else:
        edmol.RemoveBond(idx1, idx1)
    newmol = edmol.GetMol()
    #for atom in newmol.GetAtoms():
    #    if atom.GetAtomMapNum() != 0:
    #        atom.SetH
    print('remove bond',Chem.MolToSmiles(newmol))
    scaffold_smiles=''
    LG_smiles=''
    for new_smiles in Chem.MolToSmiles(newmol).split('.'):
        test_mol = Chem.MolFromSmiles(new_smiles)
        #print('test mol',new_smiles)
        for atom in test_mol.GetAtoms(): # cut and find the atom with map number is 999
            if atom.GetAtomMapNum() == 888:
                scaffold_smiles=new_smiles
            if atom.GetAtomMapNum() == 999:
                LG_smiles=new_smiles
    print('LG_smiles',LG_smiles)
    print('scaffold_smiles',scaffold_smiles)
    return LG_smiles,scaffold_smiles


def detectLeavingGroup(mol, idx, otheratom):
    if idx==-1:
        return 'H'
    atom = mol.GetAtomWithIdx(idx)
    neighbors = atom.GetNeighbors()
    if len(neighbors) == 0:
        return atom.GetSymbol()
    lv_type = atom.GetSymbol()
    for nei in neighbors:
        idx1 = nei.GetIdx()
        lv_type = lv_type + '(' + get_main(mol, idx1, idx, otherAtom)[0] + ')'
    return lv_type

data = pd.read_excel('input.xlsx')
smiles_lst=data['AAM'].values.tolist()

result = {'RXN':[], 'Type':[], 'SCF_smiles':[], 'Element':[], 'Lv_smiles':[],'subclass':[]}
sub_class_unique_list=[]
data_length=len(data)
for idx_rxn in range(data_length):
    print('row number',idx_rxn+2,'_____________________________')
    # Debug:
    #if idx_rxn < 270: continue  # <- put the id in here to debug
    #if idx_rxn == 1: break
    print(idx_rxn,data_length)
    rxn_smiles=data['AAM'][idx_rxn]
    print('rxn_smiles',rxn_smiles)
    try:
        rct_smiles,prd_smiles=rxn_smiles.split('>>')
    except:
        result['RXN'].append('')
        result['Type'].append('')
        result['SCF_smiles'].append('')
        result['Element'].append('')
        result['Lv_smiles'].append('')
        result['subclass'].append('')
        continue

    rct1s,rct2s=rct_smiles.split('.')
    rct1=Chem.MolFromSmiles(rct1s)
    rct2=Chem.MolFromSmiles(rct2s)
    rct_list=[rct1,rct2]
    prd=Chem.MolFromSmiles(prd_smiles)
    type_scf = []
    scf_smileses=[]
    res_typ = ''
    element = []
    lv_sm = []
    subclass = []
    for rct_idx in range(2):
        rct = rct_list[rct_idx]
        smiles = Chem.MolToSmiles(rct)
        print('rct',smiles)
        #print(len(data),idx_rxn)
        break_bond = data['BREAK'][idx_rxn].split('.')[rct_idx]
        print('break_bond',break_bond)
        for bond in '-=#':
            if bond in break_bond:
                break
        at_info1 = break_bond.split(bond)[0]
        at_info2 = break_bond.split(bond)[1]
        elem1 = Get_elem(at_info1)
        elem2 = Get_elem(at_info2)
        emap1 = Get_emap(at_info1)
        emap2 = Get_emap(at_info2)
        idx1 = -1
        idx2 = -1
        elem = ''
        emap = 0
        lv_elem = ''
        if emap1 == 0:
            lv_elem = elem1
            elem = elem2
            emap = emap2
        if emap2 == 0:
            lv_elem = elem2
            elem = elem1
            emap = emap1
        print('lv_elem',lv_elem)
        for atom in rct.GetAtoms():
            if atom.GetSymbol() == elem and atom.GetAtomMapNum() == emap:
                idx1 = atom.GetIdx()
                break
        print('atom',atom.GetSymbol())
        for bond in rct.GetAtomWithIdx(idx1).GetBonds():
            otherAtom = rdchem.Bond.GetOtherAtom(bond, rct.GetAtomWithIdx(idx1))
            print('*********other atom',otherAtom.GetSymbol())
            if otherAtom.GetAtomMapNum() == 0:
                idx2 = otherAtom.GetIdx()
                break
        print('atom index 1:',idx1 ,' (atom in scaffold side) map no:',atom.GetAtomMapNum())
        print('atom index 2:', idx2,'(atom in leaving group side) map no:',otherAtom.GetAtomMapNum())
        element.append(lv_elem)
        typ,scf_SMILES = get_main(rct, idx1, idx2, otherAtom) # run through each cutting atom
        type_scf.append(typ)
        scf_smileses.append(get_lg_n_scf_smiles(rct,idx1,idx2)[1])
        print('scaffold',typ)
        lv_sm.append(get_lg_n_scf_smiles(rct,idx1,idx2)[0])
        subclass.append(detectLeavingGroup(rct, idx2, otherAtom))

    if len(type_scf) == 2 :
        res_typ = type_scf[0] + '.' + type_scf[1]
        result['RXN'].append(rxn_smiles)
        result['Type'].append(res_typ)
        result['SCF_smiles'].append(scf_smileses[0]+'.'+scf_smileses[1])
        result['Element'].append(element[0] + '.' + element[1])
        result['Lv_smiles'].append(lv_sm[0] + '.' + lv_sm[1])
        result['subclass'].append(subclass[0] + '.' + subclass[1])
    #print(res_typ)
    #break #<- debug

type_list=result['Type']


# Generalize aryl:
type_list_gen=[]
for type_string in type_list:
    type_string=re.sub(r'Ar+\d', 'Ar', type_string)
    #type_string=re.sub('ArHet', 'Ar', type_string) #<- this one is wrong
    type_list_gen.append(type_string)

type_list_set=[]
type_list_out=[]
type_list_display=[]
type_smiles_list=[]
EG_class_lst=[]
EG_class_display=[]
EG_smiles_lst=[]
EG_sub_class_lst=[]

# Unique in order
arranged_smiles_lst_out=[]
for i,itm in enumerate(type_list_gen):
    # prepare smiles
    smiles=smiles_lst[i]
    rct_smiles,prd_smiles=smiles.split('>>')
    rct1_smiles,rct2_smiles=smiles.split('.')

    # prepare scaffold
    SCF_smiles=result['SCF_smiles'][i]
    #SCF_gen=result['SCF_gen'][i]
    SCFS1,SCFS2=SCF_smiles.split('.')
    pair_smiles=SCFS1+'.'+SCFS2
    pair_rev_smiles=SCFS2+'.'+SCFS1

    SCF1,SCF2=itm.split('.')
    pair=SCF1+'.'+SCF2
    pair_rev=SCF2+'.'+SCF1
    print('original pair:',pair)
    #Prepare other swapping:
    EG_class=result['Element'][i]
    EG_class=re.sub(r"H\d+","",EG_class)
    EG_class_rev=string_swap(EG_class)
    EG_sub_class=result['subclass'][i]
    EG_sub_class_rev=string_swap(EG_sub_class)
    EG_smiles=result['Lv_smiles'][i]
    EG_smiles_string=[]

    # Handle the H case
    EG_class_entry_lst=EG_class.split('.')
    EG_smiles_entry_lst=EG_smiles.split('.')
    print('EG_class_entry_lst',EG_class_entry_lst)
    print('EG_smiles_entry_lst',EG_smiles_entry_lst)
    EG_smiles_entry_lst_out=[]
    for index, item in enumerate(EG_smiles_entry_lst):
        print('index',index)
        if EG_class_entry_lst[index] == 'H':
            EG_smiles_entry_lst_out.append('[H]')
        else:
            if item=='':
                EG_smiles_entry_lst_out.append('[H]')
            else:
                EG_smiles_entry_lst_out.append(item)
    EG_smiles='.'.join(EG_smiles_entry_lst_out)
    EG_smiles=EG_smiles.replace('[[H]]','[H]')
    EG_smiles_rev=string_swap(EG_smiles)
    EG1,EG2=EG_class.split('.')
    print(EG1,EG2)
    if pair and pair_rev not in type_list_set:
        type_list_set.append(pair)
    if pair_rev in type_list_set:
        EG_class_display.append(EG1+'.'+EG2)
        type_list_display.append(pair_rev)
        #   print('unique pair:',pair_rev)
    else:
        EG_class_display.append(EG2+'.'+EG1)
        type_list_display.append(pair)
        #print('unique pair:',pair)

    EG1EN=get_en(str(EG1))
    print(EG1EN)
    EG2EN=get_en(str(EG2))
    if EG1EN>=EG2EN:
        EG_smiles_lst.append(EG_smiles)
        EG_class_lst.append(EG_class)
        EG_sub_class_lst.append(EG_sub_class)
        type_list_out.append(pair)
        type_smiles_list.append(pair_smiles)
    else:
        EG_smiles_lst.append(EG_smiles_rev)
        EG_class_lst.append(EG_class_rev)
        EG_sub_class_lst.append(EG_sub_class_rev)
        type_list_out.append(pair_rev)
        type_smiles_list.append(pair_rev_smiles)


data['SCF']=result['Type']
data['SCF_gen']=type_list_gen
data['SCF_gen_display']=type_list_display
data['SCF_smiles']=type_smiles_list
data['EG_class']=EG_class_lst
data['EG_class_display']=EG_class_display
data['EG_sub_class']=EG_sub_class_lst
data['EG_smiles']=EG_smiles_lst

df = pd.DataFrame(result)
data.to_excel('result.xlsx', index=False)
