# Load all connections as reactant -> product
in_file = open("reaction_traverse.txt",'r')
lib = {}
for line in in_file:
    line = line[:-1]
    [rea, pro] = line.split(':')
    lib[rea] = pro.split(',')
    while '' in lib[rea]:
        lib[rea].remove('')
in_file.close()

# Choose core, set up
core = ['16','21491','8','39089','256','9',
        '52','17841','317535','9676'] #top10
MAX = 4867375
visited = [False for i in range(MAX)]
queue = [x for x in core]
path = {}
for x in range(MAX):
    path[str(x)] = []
for x in core:
    path[x] = [x]

# Main algorithm
while len(queue) > 0:
    current = queue.pop(0)
    if visited[int(current)]:
        continue
    visited[int(current)] = True
    path[current].append(current)
    
    for x in lib[current]:
        if not visited[int(x)]:
            queue.append(x)
            path[x] = [sb for sb in path[current]]

# Print path from core to each chemical
out = open("pathway_from_core.csv",'w')
out.write("ID,Distance,Path\n")
for x in range(MAX):
    out.write(str(x)+','+str(len(path[str(x)])-1)+', '
              +' >> '.join(path[str(x)])+'\n')

out.close()
        

