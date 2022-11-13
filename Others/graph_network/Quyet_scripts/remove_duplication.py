file = open("reactions_numerical_form_mayduplicated.txt",'r')
all_reactions = []
for line in file:
    line = line[:-1]
    [reactant, product] = line.split(">>")
    rea = set(sorted(reactant.split(".")))
    pro = set(sorted(product.split(".")))
    all_reactions.append('.'.join(rea)+'>>'+'.'.join(pro))
file.close()

all_reactions = set(all_reactions)
out = open("reactions_numerical_form.txt",'w')
for x in all_reactions:
    out.write(x+'\n')
out.close()

