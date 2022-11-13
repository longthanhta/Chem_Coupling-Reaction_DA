import camelot
import pandas as pd
tables = camelot.read_pdf('/media/gordon/Backup Plus/Documents/Research-Prof.Su/DataMiningProject/EsayOCR/Table_Extractor/Camelot/Test1.pdf',pages="all")
print(tables)
print(len(tables))
df1 = tables[0].df
print(tables[0].parsing_report)
print(df1)
df1.to_csv("table.csv")