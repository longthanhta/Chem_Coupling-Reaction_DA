file = open("extracted_SMILES_reaxys_7M.csv",'r')
output = open("smiles_only.txt",'w')
for line in file:
    temp_list = line.split(',')
    for x in temp_list:
        if '>>' in x:
            output.write(x+'\n')
