reaction_file = open("reactions_numerical_form.txt",'r')


NUM_OF_CHEMICALS = 4867375
lib = {}
for i in range(NUM_OF_CHEMICALS):
    lib[str(i)] = set()

for line in reaction_file:
    line = line[:-1]
    [reactant, product] = line.split(">>")
    rea = reactant.split(".")
    pro = product.split(".")
    for x in pro:
        for y in rea:
            lib[x].add(y)

reaction_file.close()
out = open("reaction_reverse.txt",'w')

for i in range(NUM_OF_CHEMICALS):
    out.write(str(i) +':'+','.join(sorted(
            list(lib[str(i)]),key = len))+'\n')

out.close()

'''
out_reverse = open("reaction_reverse.txt",'w')
reverse = {}
for i in range(4867376):
    reverse[str(i)] = set()

for line in reaction_file:
    line = line[:-1]
    [reactant, product] = line.split(">>")
    rea = reactant.split(".")
    pro = product.split(".")
    for x in rea:
        for y in pro:
            reverse[x].add(y)
for i in range(1000):
    out_reverse.write(str(i) + ':\n')
for i in range(1000,4867376):
    out_reverse.write(str(i) +':'+','.join(sorted(list(reverse[str(i)]),key = len))+'\n')
out_reverse.close()
'''


    
