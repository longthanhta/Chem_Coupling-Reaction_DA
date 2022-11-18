from glob import glob

import os

import pandas as pd

import subprocess

from tqdm import tqdm

# text mining

import pytesseract

import cv2



def graph_mining(file_path):
	try:
		img = cv2.imread(file_path, cv2.IMREAD_GRAYSCALE)
		OCR_text = pytesseract.image_to_string(img)
		return OCR_text
	except:
	#	return
		return ''









directory='segmented_figs'

os.chdir(directory)


sfl = glob("*/")


print(sfl)

rsl=[] #reaction smiles
psl=[] #product smiles
ctl=[] #condition top
cbl=[] #condition bottom
cal=[] #condition all
for sf in tqdm(sfl):
	print(sf)
	os.chdir(sf)

	# text mining

	r = subprocess.run(['osra', 'reactant.png'], stdout=subprocess.PIPE)
	r_smiles=str(r.stdout.decode())
	rsl.append(r_smiles)

	try:
		r_smiles=p_smiles.split('./n')
		r_smiles='.'.join(r_smiles)
	except:
		pass



	p = subprocess.run(['osra', 'product.png'], stdout=subprocess.PIPE)
	p_smiles=str(p.stdout.decode())


	try:
		p_smiles=p_smiles.split('./n')
		p_smiles='.'.join(p_smiles)
	except:
		pass


	psl.append(p_smiles)


	# graph mining
	
	cond_top=graph_mining('cond_top.png')
	cond_bot=graph_mining('cond_bot.png')
	cond_all=graph_mining('cond.png')
	ctl.append(cond_top)
	cbl.append(cond_bot)
	cal.append(cond_all)


	print(r_smiles,p_smiles,cond_top,cond_bot,cond_all)	
	os.chdir('..')

# write the output file
os.chdir('..')
out_df=pd.DataFrame({'doi and figs':sfl,'reactants smiles':rsl,'product smiles':psl,'condition (all)':cal,'condition (above)':ctl,'condition (below)':cbl})
out_df.to_excel('final_result.xlsx')

