import tabula
 
file = "/media/gordon/Backup Plus/Documents/Research-Prof.Su/DataMiningProject/EsayOCR/Table_Extractor/tabula/Test1.pdf"
df = tabula.read_pdf(file, pages="all")
#print(df)
#print()
print(len(df))
print()
count = 0
for i in df:
    print(i)
    i.to_csv('Result_'+str(count)+'.csv',index=False)
    count = count + 1