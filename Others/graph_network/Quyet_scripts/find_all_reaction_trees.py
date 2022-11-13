reaction_file = open("reaction_reverses.txt",'r')
tree_file = open("reaction_tree.txt",'w')

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
'''

# Step 1: Focus on only one main thread, instead of having branches
lib = {}
for line in reaction_file:
    line = line[:-1]
    [pro, rea] = line.split(':')
    lib[pro] = rea.split(',')

full_trees = []
for id_ in range(4867375,999,-1):
    id_ = str(id_)
    #if len(lib[id_]) == 0:                              # ignore chemicals having no descendants
    #    continue
    can_find = False
    for tree_id_ in range(len(full_trees)):
        if full_trees[tree_id_][-1] == id_:         # if we can find the chemical at the end of some series
            can_find = True
            for child_id_ in range(len(lib[id_])):
                temp = full_trees[tree_id_]
                if child_id_ == 0:                       # append to the current tree
                    full_trees[tree_id_].append(lib[id_][child_id_])
                else:                                           # duplicate the original tree and add other children
                    added = [x for x in temp]
                    added.append(lib[id_][child_id_])
                    full_trees.append(added)
            
    if not can_find:
        full_trees.append([id_])

# store full trees
for tree in full_trees:
    tree_file.write('>>'.join(tree) + '\n')

