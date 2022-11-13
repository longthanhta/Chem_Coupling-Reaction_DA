in_file = open("reaction_traverse.txt",'r')
lib = {}
for line in in_file:
    line = line[:-1]
    [rea, pro] = line.split(':')
    lib[rea] = pro.split(',')
    while '' in lib[rea]:
        lib[rea].remove('')

in_file.close()

core = ['16']
MAX = 4867375
visited = [False for i in range(MAX)]
distance = [-1 for i in range(MAX)]
queue = [x for x in core]

while len(queue) > 0:
    current = queue.pop(0)
    if visited[int(current)]:
        continue
    visited[int(current)] = True
    
    for x in lib[current]:
        if not visited[int(x)] and x not in queue:
            queue.append(x)
            distance[int(x)] = distance[int(current)]+1
            

out = open("distance_from_core.txt",'w')
out.write("ID,Distance\n")
for x in range(MAX):
    out.write(str(x)+','+str(distance[x])+'\n')

out.close()
        

