from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.io import read
import numpy as np




#INPUT_FILES_______________________________________________________________________________________________________

# file path for atom radii empirical data
atom_radii_path = "__data_empirical/atom_radii.txt"
# file path for a _tmp_file storage
_tmp_file_loc = 'mol_files/_tmp.mol'
atom_radii_dict = {}

with open(atom_radii_path,"r") as input_:
    content = input_.read()
    for line in content.split("\n")[1:]:
        atom_radii_dict[line.split(" ")[0]] = line.split(" ")[1]




#FUNCTIONS_______________________________________________________________________________________________________
# get distance of atom from center
def get_DistFromCent(atom_cord,atom_cent_cord):
    r_2 = 0
    for i in range(len(atom_cord)):
        r_2 += (atom_cord[i]-atom_cent_cord[i])**2
    return r_2*(1/2)

# get ster3DFeat from smiles, 1st mode: simple volume; 2nd mode: filled volume
def get_Smil2Ster3DFeat(smil, mode = "Simple_Volume"):
    mol = Chem.MolFromSmiles(smil)
    # add H to mol before calculation is turned off
    #mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol,randomSeed=0xf00d) 
    #mol = Chem.RemoveHs(mol)

    # a temp file is generated for later reading
    _tmp_file=open(_tmp_file_loc,'w+')
    _tmp_file.write(Chem.MolToMolBlock(mol))
    _tmp_file.close()

    # use read in ase to read temp file as mol file
    mol = read(_tmp_file_loc)

    # get coordinate of center of atoms AND number of atom
    atom_cent_cord = np.array([0.0 ,0.0 ,0.0 ])
    atom_count = 0
    for i in range(len(mol.get_positions())):
        atom_cord = mol.get_positions()[i]
        atom_cent_cord += atom_cord
        atom_count += 1
    atom_cent_cord = atom_cent_cord/atom_count

    # For option "Simple_Volume": get maximum distance of atoms towards coordinate of center
    r_max = 0
    atom_count_dict = {}
    for i in range(len(mol.get_positions())):
        atom_cord = mol.get_positions()[i]
        elem = mol.symbols[i]
        if elem not in atom_count_dict.keys():
            atom_count_dict[elem] = 1
        else:
            atom_count_dict[elem] += 1
        r = get_DistFromCent(atom_cord,atom_cent_cord)
        if r > r_max:
            r_max = r

    # For option "Filled_Volume": cal filled volume using atom volume
    filled_volume = 0
    for elem, _count in atom_count_dict.items():
        radius = float(atom_radii_dict[elem])
        if radius == -1:
            print("ERROR, radius info not available for atom", elem)
        filled_volume += _count*radius**3
    percent_filled_vol = filled_volume

    # return cubic of maximum radius
    if mode == "Simple_Volume":
        return  r_max**3
        
    elif mode == "Filled_Volume":
        return  filled_volume



#EXAMPLES_FOR_CHECKING____________________________________________________________________________________________
#smil = '[CH2:6][CH2:5][C:3]([O:2][CH3:1])=[O:4]'
#print(get_Smil2Ster3DFeat(smil))