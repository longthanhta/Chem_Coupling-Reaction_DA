from rdkit import Chem
from rdkit.Chem import AllChem
from CORE import *
from get_Info2Feat import *     
import copy
import json




#INPUT_FILES_______________________________________________________________________________________________________

monitor_path = 'Lvg/_tmp/MOL_normal_tmp_display.o.png'
monitor_path2 = 'Lvg/_tmp/MOL_masked_tmp_display.o.png'


# Dataset for Hammett
hammett_file_path = "__data_empirical/Hammett.txt"
hammett_elec_effect_dict = {}
with open(hammett_file_path,"r") as input_:
    for line in input_.read().split("\n")[1:]:
        smil, m_hammett, p_hammett = line.split(" ")
        m_hammett = float(m_hammett) if m_hammett != "None" else None
        p_hammett = float(p_hammett) if p_hammett != "None" else None
        hammett_elec_effect_dict[smil] = {'m':m_hammett,'p':p_hammett}

hammett_sub_group_dict_path = "__data_manual/Hammett_sub_group_dict.json"
with open(hammett_sub_group_dict_path, 'r') as input_:
    hammett_sub_group_dict = json.loads(input_.read())
    hammett_sub_group_dict = {int(k):v for k,v in hammett_sub_group_dict.items()}

hammett_strict_sub_group_dict_path = "__data_manual/Hammett_strict_sub_group_dict.json"
with open(hammett_strict_sub_group_dict_path, 'r') as input_:
    hammett_strict_sub_group_dict = json.loads(input_.read())
    hammett_strict_sub_group_dict = {int(k):v for k,v in hammett_strict_sub_group_dict.items()}




# Dataset for LvgENu
lvgenu_file_path = "__data_manual/LvgENu.txt"
lvgenu_dict = {}
with open(lvgenu_file_path,"r") as input_:
    for line in input_.read().split("\n")[1:]:
        smil, nuct_elct_label = line.split(" ")
        lvgenu_dict[smil] = nuct_elct_label

lvgenu_sub_group_dict_path = "__data_manual/LvgENu_sub_group_dict.json"
with open(lvgenu_sub_group_dict_path, 'r') as input_:
    lvgenu_sub_group_dict = json.loads(input_.read())
    lvgenu_sub_group_dict = {int(k):v for k,v in lvgenu_sub_group_dict.items()}

lvgenu_strict_sub_group_dict_path = "__data_manual/LvgENu_strict_sub_group_dict.json"
with open(lvgenu_strict_sub_group_dict_path, 'r') as input_:
    lvgenu_strict_sub_group_dict = json.loads(input_.read())
    lvgenu_strict_sub_group_dict = {int(k):v for k,v in lvgenu_strict_sub_group_dict.items()}




# metal_atom list
metal_atom_list_path = "__data_empirical/metal_atom.txt"
with open(metal_atom_list_path, 'r') as input_:
    metal_atom_list = input_.read().split("\n")




#FUNCTIONS_______________________________________________________________________________________________________

# return a mol with atoms replaced by a unusual atom (118) when the scan_displ is out of range
#CHANGETO get_Mol2MaskedMol
def Mol2MaskedMol(mol,target_atom_id,scan_displ,EXCLUDE_adj_multibond = True,masking_atom_num = 118):
    # get the displacement list of mol towards target_atom_id in terms of No. of bond
    displ_list = AllChem.GetDistanceMatrix(mol)[target_atom_id]

    mol_masked = copy.deepcopy(mol)
    for a_id in range(len(displ_list)):
        # scan_displ = 0: consider atom connecting to * only (also consider C=O e.g.),
        #            = 1: consider also atoms 1 bond away from that atom
        
        # mask all bond outside scan_displ range except for multiple bond order
        if EXCLUDE_adj_multibond and displ_list[a_id] == 1+scan_displ:
            for bond in mol_masked.GetAtomWithIdx(a_id).GetBonds():
                if bond.GetBondTypeAsDouble() < 2 :
                    mol_masked.GetAtomWithIdx(a_id).SetAtomicNum(masking_atom_num)
        # mask all bond outside scan_displ range
        elif displ_list[a_id] > scan_displ:
            mol_masked.GetAtomWithIdx(a_id).SetAtomicNum(masking_atom_num)
    return mol_masked

# return a list of one hot (or the first matched functional group if return_match = True)
#CHANGETO get_Mol2Fg
def CheckFuncGroup(mol,fg_list, return_match = False):
    fg_boolean_list = []
    for i in range(len(fg_list)):
        # find a match for functional group in given mol
        fg = Chem.MolFromSmarts(fg_list[i])
        matches_copy = mol.GetSubstructMatches(fg)
        if matches_copy==():
            fg_boolean_list.append(0)
        else:
            if return_match == True:
                return fg_list[i]
            fg_boolean_list.append(1)
    
    if return_match == True:
        return None
    return fg_boolean_list


def get_SmilANDDict2Feat(   smil,   \
                            label,   \
                            output_effect_dict, \
                            sub_group_dict, \
                            strict_sub_group_dict,  \
                            other_input = None,   \
                            label_type = "element", \
                            H_default_value = 0, \
                            EXCLUDE_adj_multibond = True,   \
                            MAX_scan_displ = None):
    output = None
    if MAX_scan_displ == None:
        MAX_scan_displ = int(list(sub_group_dict.keys())[-1])

    mol = Chem.MolFromSmiles(smil)
    for target_atom in mol.GetAtoms():
        elem = target_atom.GetSymbol()
        amap = target_atom.GetAtomMapNum()

        sub_atom = None

        if label_type == "element":
            if elem == "*":
                # there should be only 1 break pt connected by only 1 atom connected
                for atom in target_atom.GetNeighbors():
                    sub_atom = atom
                    sub_atom_elem = sub_atom.GetSymbol()
                    
                    # if sub is -H, output = default
                    if sub_atom_elem == 'H':
                        output = H_default_value
                        break
                    
                    sub_atom_id = sub_atom.GetIdx()
                    # get the maximum displacement from the sub_atom
                    MAX_displ= max(AllChem.GetDistanceMatrix(mol)[sub_atom_id][1:])
                    
                    break

        elif label_type == "amap":
            if amap == 999:
                sub_atom = target_atom
                sub_atom_id = sub_atom.GetIdx()
                MAX_displ = max(AllChem.GetDistanceMatrix(mol)[sub_atom_id])

        if sub_atom != None:


            # start scanning
            for scan_displ in range(0,MAX_scan_displ+1):
                # generate a masked mol and set options for scanning
                mol_masked = Mol2MaskedMol(mol,sub_atom_id,scan_displ,EXCLUDE_adj_multibond)

                # get corresponding sub_group to scan from manual data
                sub_group_list = sub_group_dict[scan_displ].copy()
                if MAX_displ== scan_displ:
                    sub_group_list += strict_sub_group_dict[scan_displ]

                # get match group
                sub_group_match = CheckFuncGroup(mol_masked,sub_group_list,return_match = True)

                # get output if has a match, then break from loop
                if sub_group_match != None:
                    # if there is other input
                    if other_input != None:
                        output = output_effect_dict[sub_group_match][other_input]
                    else:
                        output = output_effect_dict[sub_group_match]
                    break
            
            # treat as urgent problem and terminate if output not exist
            if output == None:
                # print unmasked version
                d2d = Draw.MolDraw2DCairo(800,300)
                d2d.DrawMolecule(mol)
                png_normal = d2d.GetDrawingText()
                open(monitor_path,'wb+').write(png_normal)     
                
                # print masked version
                try:
                    d2d = Draw.MolDraw2DCairo(800,300)
                    d2d.DrawMolecule(mol_masked)
                    png_masked = d2d.GetDrawingText()
                    open(monitor_path2,'wb+').write(png_masked)     
                except:
                    d2d = Draw.MolDraw2DCairo(800,300)
                    d2d.DrawMolecule(Mol2MaskedMol(mol,sub_atom_id,scan_displ,EXCLUDE_adj_multibond,masking_atom_num = 32))
                    png_masked = d2d.GetDrawingText()
                    open(monitor_path2,'wb+').write(png_masked)     
                print("Please Check _tmp, smil=",smil,"at path",monitor_path2)
                exit()
    return output


# return the Hammett constant
#CHANGETO get_smil2ArElec_Hammett
def smiles2ArElectronic_Hammett(smil,pos):
    # change pos to p if o
    if pos == 'o':
        pos = 'p'
    elec_effect = get_SmilANDDict2Feat( smil = smil,
                                        label = "*",
                                        output_effect_dict = hammett_elec_effect_dict,
                                        sub_group_dict = hammett_sub_group_dict,
                                        strict_sub_group_dict = hammett_strict_sub_group_dict,
                                        other_input = pos,
                                        label_type = "element", 
                                        EXCLUDE_adj_multibond = True,
                                        H_default_value = 0)
    return elec_effect


# return the Hammett constant
#CHANGETO get_smil2ArElec_Hammett
def LvgSmil2NuctElctType(smil):
    nuct_elct_type = get_SmilANDDict2Feat(  smil = smil,
                                            label = 999,
                                            output_effect_dict = lvgenu_dict,
                                            sub_group_dict = lvgenu_sub_group_dict,
                                            strict_sub_group_dict = lvgenu_strict_sub_group_dict,
                                            other_input = None,
                                            label_type = "amap", 
                                            EXCLUDE_adj_multibond = True,
                                            H_default_value = 0,
                                            MAX_scan_displ = 4)
    return nuct_elct_type

#EXAMPLES_FOR_CHECKING____________________________________________________________________________________________
