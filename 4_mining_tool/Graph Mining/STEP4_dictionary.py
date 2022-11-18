import pandas as pd


data_df=pd.read_excel('final_result.xlsx')
dict_df=pd.read_excel('dict.xlsx')
ligand_list=dict_df['ligand'].values.tolist()
base_list=dict_df['base'].values.tolist()
solvent_list=dict_df['solvent'].values.tolist()
reagent_list=dict_df['reagent'].values.tolist()


lig_out=[]
sol_out=[]
base_out=[]
reg_out=[]

for index,row in data_df.iterrows():

	cond=row['condition (all)']
	cond=cond.replace('\n',' ')
	cond=cond.split(' ')


	# reset the val
	ligand=''
	solvent=''
	base=''
	reg=''
	try:
		for item in cond:
			continue
	except:
		pass

	if '' in cond:
		cond.remove('')

	cond=''.join(cond)
	print(cond)
	#for item in cond:
	#	print(item)
	#	if item in ligand_list:
	#		ligand=item
	#	if item in solvent_list:
	#		solvent=item
	#	if item in base_list:
	#		base=item
	#	if item in reagent_list:
	#		reg=item

	for item in ligand_list:
		try:
			if item in cond:
				ligand=item
		except:
			pass

	for item in solvent_list:
		try:
			if item in cond:
				solvent=item
		except:
			pass

	for item in base_list:
		try:
			if item in cond:
				base=item
		except:
			pass

	for item in reagent_list:
		try:
			if item in cond:
				reg=item
		except:
			pass

	lig_out.append(ligand)
	sol_out.append(solvent)
	base_out.append(base)
	reg_out.append(reg)


data_df['reagent']=reg_out
data_df['ligand']=lig_out
data_df['solvent']=sol_out
data_df['base']=base_out
data_df.to_excel('final_result_w_cond.xlsx')
