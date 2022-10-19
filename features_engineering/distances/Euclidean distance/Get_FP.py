from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import plotly.express as px

from threading import Thread
from multiprocessing import Process, Manager
from math import sqrt, floor

def TanimotoSimilarity(fp1, fp2,TAKEN_COLUMNS):
	#print('fp1',fp1)
	#print('fp2',fp2)
	if(len(fp1) != len(fp2)):
		print("Finger prints have different length")
		return False
	intersection = [fp1[i]*fp2[i] for i in range(len(fp1))]
	sum_intersection = sum(intersection)
	sum_union = sum(fp1)+sum(fp2)-sum_intersection
	if sum_union == 0: 
		# that means all entries of fp1 and fp2 are 0
		return 1
	else:
		tanimoto = sum_intersection/sum_union
		return tanimoto

def TanimotoDistanceForLargeFingerPrint(fp1, fp2,TAKEN_COLUMNS,RANGE_INDEX,WEIGHT):
	dist_components = []
	for i in range(len(TAKEN_COLUMNS)):
		col = TAKEN_COLUMNS[i]
		# index range of fp of column 'col'
		lower = RANGE_INDEX[i]
		upper = RANGE_INDEX[i+1]
		# distance for 'col'
		similarities=TanimotoSimilarity(fp1[lower:upper],fp2[lower:upper],TAKEN_COLUMNS)
		#print('tanimoto',similarities)
		d = 1-similarities
		dist_components.append(d*d*WEIGHT[col])

	# RMS-Tanimoto distance
	return sqrt(sum(dist_components))


def computeTanimotoDistance(indices,fps,dist_matrix,TAKEN_COLUMNS,RANGE_INDEX,WEIGHT):
	nfps = len(fps)
	for i in indices:
		print(i+1,'/',nfps)
		for j in range(i+1):
			#print(i,j)
			distance=TanimotoDistanceForLargeFingerPrint(fps[i],fps[j],TAKEN_COLUMNS,RANGE_INDEX,WEIGHT)
			dist_matrix[i*nfps+j] = distance
			dist_matrix[j*nfps+i] = dist_matrix[i*nfps+j]

if __name__ == "__main__":

	# python plot_html_RMS-Tanimoto_distance_multiprocessing.py
	# Data file, change filename and sheetname properly before running
	df = pd.read_excel('Set.xlsx')
	df = df.fillna('NA')
	#df = pd.read_excel(xls, 'Unique')	# .fillna('')

	big_lib = {}
	writer = pd.ExcelWriter('library.xlsx')
	whole_column_list=[]
	for column_name in ['GEN_CAT_SMILES','LEAVE','FORM']:
		#print(column_name)
		column = df[column_name].values
		element_list = []
		for x in column:
			if str(x) != 'nan':
				element_list.extend(x.replace(',',' ').replace('.',' ').split(' '))
		element_list = sorted(list(set(element_list)))
		while '' in element_list:
			element_list.remove('')
		big_lib[column_name] = {}
		for i in range(len(element_list)):
			big_lib[column_name][element_list[i]] = i
		#print('big_lib',big_lib)
		for item in range(0,3):
			whole_column_list.extend(element_list)
		print('whole_column_list',whole_column_list)
		pd.DataFrame({
				'Element':element_list,
				'ID':list(range(len(element_list)))
			}).to_excel(writer,f'element_list_{column_name}')

	writer.save()

	#===============[Generate Fingerprint]=================#
	rxns = df
	NUM_ROW = len(rxns)
	fps = []

	FP_LENGTH = {}
	for col in ['GEN_CAT_SMILES','LEAVE','FORM']:
		FP_LENGTH[col] = len(big_lib[col])
	#print(FP_LENGTH)
	MAX_NUM_PARTS = {'GEN_CAT_SMILES':3,'LEAVE':3,'FORM':3} # maximum parts taken for each column, eg. 5 rcts in GEN_CAT_SMILES
	
	# Take keywords in the following list to build fingerprint ['ATOMS_GEN_CAT_SMILES','GEN_CAT_SMILES','LEAVE','CHANGE','FORM']
	TAKEN_COLUMNS = ['GEN_CAT_SMILES','LEAVE','FORM']
	WEIGHT = {'GEN_CAT_SMILES':1,'LEAVE':3,'FORM':9} # sum of those numbers whose column are taken must be 0
	# range of index, as fp for each taken column
	RANGE_INDEX = [0]
	for col in TAKEN_COLUMNS:
		RANGE_INDEX.append(RANGE_INDEX[-1]+FP_LENGTH[col]*MAX_NUM_PARTS[col])
	#print(RANGE_INDEX)

	imcomplete_rxns = []
	for i in range(NUM_ROW):
		print(f'Extract finger print: {i}/{NUM_ROW}')
#		if rxns['FORM_count'][i] > 2 or rxns['FORM_count'][i] == 0:
#			print('Number of formed bond no reasonable, ignore!')
#			imcomplete_rxns.append(i)
#			continue

		fp = []
		# finger print for each column
		for col in  TAKEN_COLUMNS:
			cell = str(rxns[col][i])
			elements = cell.split('.')
			elements.extend(['']*MAX_NUM_PARTS[col])
			for j in range(MAX_NUM_PARTS[col]):
				f = [0]*FP_LENGTH[col]
				l = elements[j].split(',')
				for e in l:
					#print('col',col)
					if e == '' : continue
					if e == 'nan' :continue
					f[big_lib[col][e]] = 1*WEIGHT[col] # if big_lib is python lib, remove ['ID']
				fp.extend(f)
				#print('fp for category',col,f)
		
		# add fp to list
		#print(fp)
		fps.append(fp)



	#==============[Save finger print]=====================#
	print("Check length of fingerprint:",len(fps),len(fps[0]))
	df = rxns.drop(imcomplete_rxns)
	df.index = range(len(df)) # reindex df from 1

	# fps = [[0,0] for i in range(fps)]
	#nfps = len(fps)
	#fps = fps[:nfps]
	
	### Calculate distance using multiprocessing
	#batch = 8 # = num of cores in computer
	#cutoff = [floor(nfps*sqrt(i/batch)) for i in range(batch+1)]
	#print(cutoff)
#
	### Multiprocessing
	#dist_matrix = Manager().list([1]*(nfps*nfps)) # 1-dim
	#process = []
#
	#for i in range(batch):
	#	p = Process(target=computeTanimotoDistance, 
	#		args=(range(cutoff[i],cutoff[i+1]),fps,dist_matrix,TAKEN_COLUMNS,RANGE_INDEX,WEIGHT))
	#	p.start()
	#	process.append(p)
	#for p in process:
	#	p.join()


	## convert dist_matrix to 2 dim list
	#print('dists',dists)
	print('whole_column_list len',len(whole_column_list))
	fps_np=np.array(fps)
	fps_df=pd.DataFrame(fps,columns=whole_column_list)
	fps_df.to_csv('fingerprint_check.csv')
	pickle.dump( fps_np, open("fingerprint.pkl", "wb" ))
	df.to_csv('data_for_plot.csv')
