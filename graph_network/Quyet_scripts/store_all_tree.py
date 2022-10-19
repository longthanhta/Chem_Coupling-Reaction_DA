tree_file = open("reaction_tree.txt",'w')
reaction_file = open("reaction_reverse.txt",'r')
dict_file = open("chem_dictionary.txt",'r')

lib_chem = {}
inverse = {}
for line in dict_file:
    line = line[:-1] # eliminate '\n'
    [id_ ,name_] = line.split(',')
    lib_chem[name_] = id_ # from smiles to id_
    inverse[id_] = name_    # from id_ to smiles
inverse['cyclic'] = ''

lib_react = {}
for line in reaction_file:
    line = line[:-1] # eliminate '\n'
    [pro, rea] = line.split(':')
    lib_react[pro] = rea.split(',')

used = [False for x in range(4867376)]

for x in range(4867374,4867373,-1):
    if used[x]: continue
    used[x] = True
    series = [[str(x)]]
    while True:
        end = True
        for thread in series:
            if thread[-1] != 'cyclic' and lib_react[thread[-1]][0] != '':
                end = False
                for child in lib_react[thread[-1]]:
                    added = [ele for ele in thread]
                    if child not in added:
                        added.append(child)
                        used[int(child)] = True
                    else:
                        added.append('cyclic')
                    series.append(added)
                series.remove(thread)

        if end:
            break
    for thread in series:
        tree_file.write('\n'+'>>'.join(thread))


