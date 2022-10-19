# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

df = pd.read_csv("extracted_SMILES_reaxys.csv", dtype = str)
data = df[df.Catalysts != " "].loc[:,["Catalysts","ID"]]
data = data.drop_duplicates()
catalysts = data.Catalysts.str.split('\n')
for i in catalysts:
    i.pop()
    noNickel = []
    for j in i:
        if (j.find('Ni') == -1 and j.find('nickel') == -1):
            noNickel.append(j)
    for j in noNickel:
        i.remove(j)

newdata = np.empty([catalysts.count()], object)
cnt = 0
for i in catalysts.index:
    print('extract_catalyst',i)
    newdata[cnt] = pd.DataFrame({"Catalysts":catalysts[i],"ID":data.ID[i]})
    cnt += 1
newdata = pd.concat(newdata)
newdata.index = range(len(newdata))
newdata.Catalysts = newdata.Catalysts.str.strip()
newdata = newdata.drop_duplicates()
unique = newdata.drop_duplicates(subset = "Catalysts")
counts = pd.Series(np.zeros(len(unique)), index = unique.index ,name = "counts")
             
for i in unique.index:
    print('counting unique and removing duplicates',i)
    list1 = [unique.loc[i,'ID']]
    for j in newdata.index:
        if unique.Catalysts[i] == newdata.Catalysts[j]:
            counts[i] += 1
            str1 = newdata.loc[j,'ID']
            if (str1 not in list1):
                list1.append(str1)
    unique.loc[i,'ID'] = pd.Series(list1).str.cat(sep = ',')

final = pd.DataFrame({"Catalysts":unique.Catalysts,"ID":unique.ID,"counts":counts})
final.index = range(len(final))
##print(final)
final.to_csv("extract_catalyst_2.csv")