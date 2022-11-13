with open('SMARTS_template.txt', 'r') as f:
    data_temp = f.readlines()
template_ID = input("input template ID (not including the first numbe 1) ")
template_in=data_temp[int(template_ID)]
print('SMARTS',template_in)

import pandas as pd
import json
import ast
with open('RAW.txt', 'r') as f:
    data = f.readlines()
template=str(template_in).rstrip()
length=len(data)
for i in range(0,length):
	if i%10000==0:
		print('searching progress',str(i*100/length),'%')
	line_string=str(data[i])
	dict_line = ast.literal_eval(line_string)
	template_file=dict_line['reaction_smarts'].rstrip()
	if template == template_file:
		ref_list=dict_line['references']
		data_react=pd.read_csv('SMILES_reaction.csv')
		print(data_react['ID'])
		print('first in ref list',ref_list[0])
		founded_ID_list=[]
		for ID in ref_list:
			ID,number=ID.split('-')
			if ID in founded_ID_list: continue
			founded_ID_list.append(ID)
			ID=int(ID)
			data_react=data_react[(data_react['ID'] == ID)]
			if len(data_react)!=0:
				print('ref for reaction',data_react['reaction'].tolist())
				DOI_list=data_react['DOI'].tolist()
				for item in DOI_list:
					print(item)
		break
