file = open("smiles_only.txt",'r')
dict_file = open("chem_dictionary.txt",'w')

set_chemicals = set()
for line in file:
    line = line[1:-2]
    [reactant, product] = line.split(">>")
    rea = reactant.split(".")
    pro = product.split(".")
    for chemical in rea:
        if(chemical != ''):
            set_chemicals.add(chemical)
    for chemical in pro:
        if(chemical != ''):
            set_chemicals.add(chemical)


list_chemicals = sorted(list(set_chemicals), key = len)
for i in range(len(list_chemicals)):
    dict_file.write(str(i) + ',' + list_chemicals[i] + '\n')

dict_file.close()
file.close()
