import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Draw
import csv

# input
n_limit = 500
task_description = "_rct2_500n" #rct = reactant; 2 = two are involved; 500n = 500 as maximum node
input_type = "reactant"
image_output_path = "images"
excel_input_path = '../reaction_fixed_precatalyst_problem.xlsx'

# read excel file
df = pd.read_excel(excel_input_path)

# dataframe for only reactant
df_input = df[['RCT','PRD']]

# convert to a list of tuple
list_input = df_input.reset_index().values.tolist()
rct_list = [i[1] for i in list_input]
prd_list = [i[2] for i in list_input]

if input_type == "reactant":
    inputLabel_list = rct_list
elif input_type == "product":
    inputLabel_list = prd_list

# function listed for clarification
def SMILES2FP(SMILES,image_path = None, index=None):
    mol = Chem.MolFromSmiles(SMILES)
    if image_path and index:
        Draw.MolToFile(mol,image_path+'/mol{}.png'.format(index)) 
    return AllChem.GetMorganFingerprintAsBitVect(mol,2,invariants=[1]*mol.GetNumAtoms())

def Tanimoto_Sim(FP_1,FP_2):
    return DataStructs.FingerprintSimilarity(FP_1,FP_2, metric=DataStructs.DiceSimilarity)

# generate input fingerprint list
inputFP_list = []
for i in range(len(inputLabel_list)):
    inputFP_list.append(SMILES2FP(inputLabel_list[i],"images",i))
    if i>n_limit-2:
        break
print("Total length of FP entries =",len(inputFP_list))

edge_list = []
weights = {}

if n_limit == None:
    n_limit = len(inputFP_list)
else:
    inputLabel_list = inputLabel_list[:n_limit]

import time
start = time.time()
# weight and edge_list generator (non-directed graph)
for i in range(n_limit):
    for j in range(n_limit):
        # Do not add parallel edges here, to be sure
        # to have the right weight later
        if i in weights and j in weights[i] or j in weights and i in weights[j]:
            continue

        # use tanimoto similarity to generate weight
        weight = Tanimoto_Sim(inputFP_list[i],inputFP_list[j])
        edge_list.append([i, j, weight])

        # Store the weights in 2d map for easy access
        if i not in weights:
            weights[i] = {}
        if j not in weights:
            weights[j] = {}

        # Invert weights to make lower ones more visible in the plot
        weights[i][j] = 1.0 - weight
        weights[j][i] = 1.0 - weight
print( "time taken for input_file=",time.time() - start,"s")

# disabled input file with visible data
'''import json

start = time.time()
with open("tmap_input{}.txt".format(task_description),"w+") as output:
    json.dump({"n_limit":n_limit,"edge_list":edge_list,"weights":weights,"label_list":inputLabel_list}, output)
print( "time taken =",time.time() - start,"s")'''

# pickled input file
import pickle
outfile = open("tmap_input{}".format(task_description),'wb')
pickle.dump({"n_limit":n_limit,"edge_list":edge_list,"weights":weights,"label_list":inputLabel_list},outfile)
outfile.close()
