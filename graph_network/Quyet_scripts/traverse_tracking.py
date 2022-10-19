reaction_file = open("reactions_numerical_form.txt",'r')
out = open("reaction_traverse.txt",'w')

NUM_OF_CHEMICALS = 4867375
lib = {}
for i in range(NUM_OF_CHEMICALS):
    lib[str(i)] = set()

for line in reaction_file:
    line = line[:-1]
    [reactant, product] = line.split(">>")
    rea = reactant.split(".")
    pro = product.split(".")
    for x in rea:
        for y in pro:
            lib[x].add(y)

for i in range(NUM_OF_CHEMICALS):
    out.write(str(i) +':'+','.join(sorted(
            list(lib[str(i)]),key = len))+'\n')

out.close()
reaction_file.close()



    
