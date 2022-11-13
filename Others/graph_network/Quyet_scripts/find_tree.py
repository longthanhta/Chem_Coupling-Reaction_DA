reaction_file = open("reaction_reverse.txt",'r')
dict_file = open("chem_dictionary.txt",'r')

''' Format of a tree
[Pro,[React (now become a Pro)[...]]]
But before that, we try something simplier
Step 1: Simple tree, e.g 1234>>563476>>23454345>>...
    For chemicals having ID < 1000, we consider them as able to be modulated
    Chemicals with no followed chemicals in Reaction_Dictionary are terminal element
        in reaction series
    Algorithm: if we can find the current chemical at the end of some series,
                just append all of its children to those series (need duplicating),
                otherwise creat a new tree
It seem not working well. Maybe we need to modify it as: enter a chemical, translate it
    to id_, then find all thread with this chemicals as the head
    - Find a thread
    - Find all threads
    - Loop all over chemicals
    - Store all trees in a file
Unfortunately, we face memory issue. Therefore, the algorithm has to change as
find the child with the largest number of childs
'''
lib_chem = {}
inverse = {}
for line in dict_file:
    line = line[:-1] # eliminate '\n'
    [id_ ,name_] = line.split(',')
    lib_chem[name_] = id_
    inverse[id_] = name_
inverse['cyclic'] = ''

lib_react = {}
for line in reaction_file:
    line = line[:-1]
    [pro, rea] = line.split(':')
    lib_react[pro] = rea.split(',')
    
x = input("Enter the chemical: ")
series = [[lib_chem[x]]]
while True:
    end = True
    for thread in series:
        if thread[-1] != 'cyclic' and lib_react[thread[-1]][0] != '':
            end = False
            for child in lib_react[thread[-1]]:
                added = [x for x in thread]
                if child not in added:
                    added.append(child)
                else:
                    added.append('cyclic')
                series.append(added)
            series.remove(thread)

    if end:
        break
for thread in series:
    print('\n'+'>>'.join([inverse[x] for x in thread]))
    #print('\n'+'>>'.join(thread))
