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

##=============== Plot function
def plot_scatter(x,df):
    fig = px.scatter(df, x=x[:,0], y=x[:,1],
        color=df['FORM(upper)'].astype(str),hover_data=['REF','AUT','ID','BREAK','LEAVE','CHANGE','FORM','Count'])
    fig.write_html("plot_manual_clusters.html")




if __name__ == "__main__":

	# python plot_html_RMS-Tanimoto_distance_multiprocessing.py
	# Data file, change filename and sheetname properly before running
	xls = pd.ExcelFile('7_Aug_labeled_N_to_Ni_bond_extraction_ordered_CAT.xlsx')
	df = pd.read_excel(xls, 'test')	# .fillna('')

	#===========[ Add and normalize data ]===========#
	df['FORM(upper)'] = df['FORM'].values

	#for i in range(len(df)):
	#	for col in ['BREAK','CHANGE','LEAVE','FORM','FORM(upper)']:
	#		df[col][i] = normalize(df[col][i],col)

	# add ATOM column
	#atoms_break = []
	#for row in df['BREAK']:
	#	if str(row) == 'nan':
	#		atoms_break.append('')
	#	else:
	#		for c in '.,=-#:':
	#			row = row.replace(c,' ')
	#		atoms_break.append(','.join(set(row.split(' '))))
	#df['ATOMS_BREAK'] = atoms_break # atom for all rcts, not just for one

	#======[Create entity libraries and write out]======#
	big_lib = {}
	writer = pd.ExcelWriter('library.xlsx')

	for column_name in ['BREAK','LEAVE','CHANGE','FORM(upper)']:
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
	for col in ['BREAK','LEAVE','CHANGE','FORM(upper)']:
		FP_LENGTH[col] = len(big_lib[col])
	#print(FP_LENGTH)
	MAX_NUM_PARTS = {'BREAK':5,'LEAVE':5,'CHANGE':5,'FORM(upper)':5} # maximum parts taken for each column, eg. 5 rcts in BREAK
	
	# Take keywords in the following list to build fingerprint ['ATOMS_BREAK','BREAK','LEAVE','CHANGE','FORM(upper)']
	TAKEN_COLUMNS = ['BREAK','LEAVE','CHANGE','FORM(upper)']
	WEIGHT = {'BREAK':1,'LEAVE':1,'CHANGE':1,'FORM(upper)':10} # sum of those numbers whose column are taken must be 0
	# range of index, as fp for each taken column
	RANGE_INDEX = [0]
	for col in TAKEN_COLUMNS:
		RANGE_INDEX.append(RANGE_INDEX[-1]+FP_LENGTH[col]*MAX_NUM_PARTS[col])
	#print(RANGE_INDEX)

	imcomplete_rxns = []
	for i in range(NUM_ROW):
		print(f'Extract finger print: {i}/{NUM_ROW}')
		if rxns['FORM_count'][i] > 2 or rxns['FORM_count'][i] == 0:
			print('Number of formed bond no reasonable, ignore!')
			imcomplete_rxns.append(i)
			continue

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
					f[big_lib[col][e]] = 1 # if big_lib is python lib, remove ['ID']
				fp.extend(f)
				#print('fp for category',col,f)
		
		# add fp to list
		fps.append(fp)


	#==============[Save finger print]=====================#
	print("Check length of fingerprint:",len(fps[0]),len(fps[2]))
	df = rxns.drop(imcomplete_rxns)
	df.index = range(len(df)) # reindex df from 1

	''' Save fingerprint and data to plot
	df.to_csv('data_for_plot.csv')
	file = open('newfingerprint','wb')
	pickle.dump(fps,file)
	'''

	## Load fingerprint
	'''filehandler = open('newfingerprint', 'rb')
	fps = pickle.load(filehandler)
	filehandler.close()'''

	# fps = [[0,0] for i in range(fps)]
	nfps = len(fps)
	fps = fps[:nfps]
	
	## Calculate distance using multiprocessing
	batch = 8 # = num of cores in computer
	cutoff = [floor(nfps*sqrt(i/batch)) for i in range(batch+1)]
	print(cutoff)

	## Multiprocessing
	dist_matrix = Manager().list([1]*(nfps*nfps)) # 1-dim
	process = []

	for i in range(batch):
		p = Process(target=computeTanimotoDistance, 
			args=(range(cutoff[i],cutoff[i+1]),fps,dist_matrix,TAKEN_COLUMNS,RANGE_INDEX,WEIGHT))
		p.start()
		process.append(p)
	for p in process:
		p.join()


	## convert dist_matrix to 2 dim list
	dists = []
	for i in range(nfps):
		dists.append(dist_matrix[i*nfps:(i*nfps+nfps)])
	print('dists',dists)

	pickle.dump( dists, open("distance_matrix.pkl", "wb" ))
	df.to_csv('data_for_plot.csv')
