in_file = open("reaction_tree_2.txt",'r')
out_file = open("reaction_tree_SMILES.txt",'w')
dict_file = open("chem_dictionary.txt",'r')

#lib_chem = {}
inverse = {}
for line in dict_file:
    line = line[:-1]               # eliminate '\n'
    [id_ ,name_] = line.split(',')
    #lib_chem[name_] = id_ # from smiles to id_
    inverse[id_] = name_    # from id_ to smiles
inverse['cyclic'] = 'cyclic'

for line in in_file:
    line = line[:-1]
    series = line.split(">>")
    for i in range(len(series)):
        series[i] = inverse[series[i]]
    out_file.write(">>".join(series) + '\n')

in_file.close()
out_file.close()
dict_file.close()
