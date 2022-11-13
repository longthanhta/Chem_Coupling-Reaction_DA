file = open("smiles_only.txt",'r')
dict_file = open("chem_dictionary.txt",'r')
out = open("reactions_numerical_form.txt",'w')

lib = {}
for line in dict_file:
    line = line[:-1] # eliminate '\n'
    [id_ ,name_] = line.split(',')
    lib[name_] = id_

for line in file:
    line = line[1:-2]
    [reactant, product] = line.split(">>")
    rea = reactant.split(".")
    pro = product.split(".")
    temp = '.'.join(set([lib[x] for x in rea])) + ">>" + '.'.join(set([lib[x] for x in pro]))
    out.write(temp + '\n')

file.close()
dict_file.close()
out.close() 
