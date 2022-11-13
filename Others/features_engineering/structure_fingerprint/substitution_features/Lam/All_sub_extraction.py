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



def check_phenyl(mol, ac_i):
    whole_system = GetRingSystems(mol)
    print('check_phenyl_whole_system: ', whole_system)
    core_system = get_core_ring(ac_i, whole_system)
    print('check_phenyl_core_system: ', core_system)
    if len(core_system) == 6:
        return 0
    elif len(core_system) > 6:
        return 1
    return -1


def check_ring_aromatic(system, mol):
    for i in system:
        at = mol.GetAtomWithIdx(i)
        if not at.GetIsAromatic():
            return False
    return True

def GetRingSystems(mol):
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and check_ring_aromatic(ringAts, mol) and check_ring_aromatic(system, mol):
                ringAts = ringAts.union(system)
            else:
                if check_ring_aromatic(system, mol):
                    nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems
def get_two_ring(mol, core_ring):
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        if len(ringAts.intersection(core_ring)):
            systems.append(ringAts)
    return systems

def get_core_ring(ac_i, system):
    for sys in system:
        if ac_i in sys:
            return sys
    return
def get_second_ring(ac_i, system):
    for sys in system:
        if ac_i not in sys:
            return sys
    return
def get_common_atom(system):
    return system[0].intersection(system[1])

def get_core_atom(SCF_mol):
    ac=False
    ac_i=-1 # set -1 because 0 also an index
    for atom in SCF_mol.GetAtoms():
        am=atom.GetAtomMapNum()
        if am == 888: # this value is pre defined in the scaffold extraction script
            ac=atom # ac means atom core
            ac_i=atom.GetIdx()
            break
    return ac, ac_i
def get_core_atom_for_second_ring(ac, core_ring, common_atom):
    res = 0
    for b1 in ac.GetBonds():
        oa1=rdchem.Bond.GetOtherAtom(b1,ac) 
        oa1_i=oa1.GetIdx()
        if oa1_i in core_ring:
            if oa1_i in common_atom:
                return oa1_i
            for b2 in oa1.GetBonds():
                oa2=rdchem.Bond.GetOtherAtom(b2,oa1)
                oa2_i=oa2.GetIdx()
                if oa2_i in common_atom:
                    res = oa2_i
    return res


def detect_sub(at_i, oa_i, mol):
    edmol = Chem.EditableMol(mol)
    edmol.RemoveBond(at_i, oa_i)
    newmol = edmol.GetMol()
    nm_s_all=Chem.MolToSmiles(newmol) # nm : new mol
    nm_s_lst=nm_s_all.split('.')
    # remove the main ring
    print('nm_s_lst',nm_s_lst)
    for nm_s in nm_s_lst:
        if '888' in nm_s: continue
        #print('atom index',ra_i,'with position',position)
        print('has substitution',nm_s)
        nm_s = re.sub(r'\:\d+', '', nm_s) # remove atom map
        #nm_s = re.sub(r'\:\d+', '', nm_s)
        print(nm_s)
        return nm_s
    return ''
        

def detect_position(ac, ac_i, ring_system, core_ring, mol, order):
    ortho = []
    meta = []
    para = []
    for b1 in ac.GetBonds():
        oa1 = rdchem.Bond.GetOtherAtom(b1,ac)
        oa1_i = oa1.GetIdx()
        if oa1_i in core_ring:
            ortho.append(oa1_i)
    for ai2 in ortho:
        at = mol.GetAtomWithIdx(ai2)
        for b2 in at.GetBonds():
            oa2 = rdchem.Bond.GetOtherAtom(b2,at)
            oa2_i = oa2.GetIdx()
            if oa2_i in core_ring and oa2_i != ac_i:
                meta.append(oa2_i)
    for ai3 in meta:
        at = mol.GetAtomWithIdx(ai3)
        for b3 in at.GetBonds():
            oa3 = rdchem.Bond.GetOtherAtom(b3, at)
            oa3_i = oa3.GetIdx()
            if oa3_i in core_ring and oa3_i not in ortho:
                para.append(oa3_i)
    sub = []
    pos = []
    position = ''
    for at_i in ortho:
        at = mol.GetAtomWithIdx(at_i)
        for b in at.GetBonds():
            oa = rdchem.Bond.GetOtherAtom(b, at)
            oa_i = oa.GetIdx()
            if oa_i not in ring_system:
                sub.append(detect_sub(at_i, oa_i, mol))
                position = 'o'+order
                pos.append(position)
    for at_i in meta:
        at = mol.GetAtomWithIdx(at_i)
        for b in at.GetBonds():
            oa = rdchem.Bond.GetOtherAtom(b, at)
            oa_i = oa.GetIdx()
            if oa_i not in ring_system:
                sub.append(detect_sub(at_i, oa_i, mol))
                position = 'm'+order
                pos.append(position)
    check =False
    for at_i in para:
        if check:
            break
        at = mol.GetAtomWithIdx(at_i)
        for b in at.GetBonds():
            oa = rdchem.Bond.GetOtherAtom(b, at)
            oa_i = oa.GetIdx()
            if oa_i not in ring_system:
                check = True
                sub.append(detect_sub(at_i, oa_i, mol))
                postion = 'p'+order
                pos.append(postion)
                
    print('ortho: ' ,ortho)
    print('meta: ' , meta)
    print('para: ' ,para)
    print('sub: ', sub)
    print('pos: ', pos)
    return sub, pos

def process_phenyl(ac, ac_i, mol, order=''):
    ring_sys = GetRingSystems(mol)
    core_ring = get_core_ring(ac_i, ring_sys)
    sub_temp, pos_temp = detect_position(ac, ac_i, core_ring, core_ring, mol, order)
    sub = []
    pos = []
    for s in sub_temp:
        sub.append(s)
    for p in pos_temp:
        pos.append(p)
    return sub, pos
def process_naph(ac, ac_i, mol, order='s'):
    ring_sys = GetRingSystems(mol)
    core_ring_sys = get_core_ring(ac_i, ring_sys)
    two_ring = get_two_ring(mol, core_ring_sys)
    core = get_core_ring(ac_i, two_ring)
    sec = get_second_ring(ac_i, two_ring)
    sub1, pos1 = detect_position(ac, ac_i, core_ring_sys, core, mol, '')
    common_atom = get_common_atom(two_ring)
    coa2i = get_core_atom_for_second_ring(ac, core, common_atom)
    coa2 = mol.GetAtomWithIdx(coa2i)
    sub2, pos2 = detect_position(coa2, coa2i, core_ring_sys, sec, mol, order)
    sub = []
    pos = []
    for s in sub1:
        sub.append(s)
    for s in sub2:
        sub.append(s)
    for p in pos1:
        pos.append(p)
    for p in pos2:
        pos.append(p)
    return sub, pos

def process_benzylic(ac, ac_i, mol):
    cons_lst = []
    sub = []
    pos = []
    order = 1
    for b1 in ac.GetBonds():
        cons=rdchem.Bond.GetOtherAtom(b1,ac) # oa means other tom
        cons_archeck=cons.GetIsAromatic()
        cons_i=cons.GetIdx()
        if cons_archeck and cons_i not in cons_lst:
            cons_lst.append(cons_i)
            if check_phenyl(mol, cons_i) == 0:
                sub_temp, pos_temp = process_phenyl(cons, cons_i, mol, str(order))
            elif check_phenyl(mol, cons_i) == 1:
                sub_temp, pos_temp = process_naph(cons, cons_i, mol, str(order) + 's')
            else: 
                continue
            order+=1
            for temp in sub_temp:
                sub.append(temp)
            for temp in pos_temp:
                pos.append(temp)
    return sub, pos

def position_detection(SCF_smiles_input):
    # Get Atom Core

    print('SCF_smiles_input_________________________________',SCF_smiles_input)
    mol=Chem.MolFromSmiles(SCF_smiles_input)

    ac, ac_i = get_core_atom(mol)
    print('ac_i: ', ac_i)
    print('ac_map', ac.GetAtomMapNum())
    
    if ac.GetIsAromatic():
        if check_phenyl(mol, ac_i) == 0:
            sub, pos = process_phenyl(ac, ac_i, mol)
            sub_s = '.'.join(sub)
            pos_s = '.'.join(pos)
            return sub_s, pos_s
        elif check_phenyl(mol, ac_i) == 1:
            sub, pos = process_naph(ac, ac_i, mol)
            sub_s = '.'.join(sub)
            pos_s = '.'.join(pos)
            return sub_s, pos_s
    sub, pos = process_benzylic(ac, ac_i, mol)
    sub_s = '.'.join(sub)
    pos_s = '.'.join(pos)
    return sub_s, pos_s

        #print(no_bond)
if __name__ == '__main__':
    data_df=pd.read_excel('benzylic.xlsx')
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
    data_df.to_excel('outall.xlsx',index=None)