in_file = open("connections.txt",'r')
dictionary = open("chem_dictionary.txt",'r')

NUM_OF_CHEMICALS = 4867375
MAX_NUM_OF_CONNECTION = 57000
n10 = n100 = n1000 = 0
arr = []
dic = {}
num = [0 for i in range(MAX_NUM_OF_CONNECTION)]

for i in range(NUM_OF_CHEMICALS):
    
    line = dictionary.readline()
    line = line[:-1]
    [id_,smile_] = line.split(',') # id_ = str(i)
    dic[id_] = smile_
    arr.append([id_,0])
    
    line = in_file.readline()
    line = line[:-1]
    [id_,noc_] = line.split(',')
    noc_ = int(noc_)
    arr[i][1] = noc_
    num[noc_] += 1
    
    if noc_ > 10: 
        n10+=1
        if noc_ > 100: 
            n100+=1
            if noc_ > 1000: 
                n1000+=1

in_file.close()
print(n10)
print(n100)
print(n1000)
dictionary.close()
arr.sort(key = lambda pair:pair[1])

out = open("counting_by_num_of_connection.txt",'w')
top30 = arr[-30:]
for i in range(10):
    out.write(top30[i][0]+','+dic[top30[i][0]]+': '+
                  str(top30[i][1])+" connections\n")
out.close()

out = open("chem_arranged_by_numOfConnection.txt",'w')
for i in range(NUM_OF_CHEMICALS):
    out.write(arr[NUM_OF_CHEMICALS-1-i][0]+','+
              dic[arr[NUM_OF_CHEMICALS-1-i][0]]+': '+
                  str(arr[NUM_OF_CHEMICALS-1-i][1])+" connections\n")
out.close()

out = open("distribution.csv",'w')
out.write("NumOfConnection,NumOfSubstance\n")
for i in range(MAX_NUM_OF_CONNECTION):
    out.write(str(i)+","+str(num[i])+'\n')
out.close()







