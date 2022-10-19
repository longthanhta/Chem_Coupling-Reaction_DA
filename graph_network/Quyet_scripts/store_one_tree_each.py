tree_file = open("reaction_tree_2.txt",'w')
reaction_file = open("reaction_reverse.txt",'r')
dict_file = open("chem_dictionary.txt",'r')

# Load dictionary file to convert id_ to SMILE and vice versa
lib_chem = {}
inverse = {}
for line in dict_file:
    line = line[:-1] # eliminate '\n'
    [id_ ,name_] = line.split(',')
    lib_chem[name_] = id_ # from smiles to id_
    inverse[id_] = name_    # from id_ to smiles
inverse['cyclic'] = ''

# Load reverse reactions in format of dictionary: lib[product] = list of possible reactants
lib_react = {}
for line in reaction_file:
    line = line[:-1] # eliminate '\n'
    [pro, rea] = line.split(':')
    lib_react[pro] = rea.split(',')
lib_react[''] = []

used = [False for x in range(4867376)] # to mark the used chemicals, not use again

for x in range(4867374,999,-1):
    if used[x]: continue
    used[x] = True
    thread = [str(x)] # only one thread, this is the reverse reaction chain with starting node as x (id_)
    while True: # this loop step-by-step add 1 chemical to the chain
        end = True
        if thread[-1] != 'cyclic' and lib_react[thread[-1]][0] != '' and len(thread) < 20: # restrict the length, CAN MODIFY to 2 to get pair reactions as OUR TASK today!
            added_child = '';
            for child in lib_react[thread[-1]]: # choose the child with the largest number of childs
                if len(lib_react[child]) > len(lib_react[added_child]) and child not in thread:
                    added_child = child
            if added_child != '':
                thread.append(added_child)
                if not used[int(added_child)]:
                    # used[int(added_child)] = True
                    end = False         
                   
        if end:
            break
    print(x)
    tree_file.write('>>'.join(thread)+'\n')

# 2 hours to run
