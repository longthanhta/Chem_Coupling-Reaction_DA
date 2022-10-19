rea_to_pro = open("reaction_traverse.txt",'r')
pro_to_rea = open("reaction_reverse.txt",'r')
dictionary = open("chem_dictionary.txt",'r')

NUM_OF_CHEMICALS = 4867375
arr = []
dic = {}
for i in range(NUM_OF_CHEMICALS):
    
    line = dictionary.readline()
    line = line[:-1]
    [id_,smile_] = line.split(',') # id_ = str(i)
    dic[id_] = smile_
    arr.append([id_,[0,0,0]])# connections as pro, rea, total
    
    line = pro_to_rea.readline()
    line = line[:-1]
    [s,line] = line.split(':') # s = id_ 
    if line != '':
        arr[i][1][2] += line.count(',')+1
    
    line = rea_to_pro.readline()
    line = line[:-1]
    [s,line] = line.split(':') # s = id_ 
    if line != '':
        arr[i][1][1] += line.count(',')+1
    arr[i][1][0] = arr[i][1][1]+arr[i][1][2]

rea_to_pro.close()
pro_to_rea.close()

arr.sort(key = lambda pair:pair[1][0])
out_order = open("chem_NOC.csv",'w')
out_order.write(
        "ID,Smile,Total,Reactions as Pro, Reactions as Rea")
for i in range(NUM_OF_CHEMICALS):
    out_order.write(arr[i][0]+','+dic[arr[i][0]]+','+
                  str(arr[i][1][0])+','+str(arr[i][1][2])+
                  ','+str(arr[i][1][1])+"\n")

out_order.close()
dictionary.close()









