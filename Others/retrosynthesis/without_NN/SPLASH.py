#	 
#      ___                         ___       ___           ___           ___     
#     /  /\          ___          /  /\     /  /\         /  /\         /  /\    
#    /  /::\        /  /\        /  /:/    /  /::\       /  /::\       /  /:/    
#   /__/:/\:\      /  /::\      /  /:/    /  /:/\:\     /__/:/\:\     /  /:/     
#  _\_ \:\ \:\    /  /:/\:\    /  /:/    /  /::\ \:\   _\_ \:\ \:\   /  /::\ ___ 
# /__/\ \:\ \:\  /  /::\ \:\  /__/:/    /__/:/\:\_\:\ /__/\ \:\ \:\ /__/:/\:\  /\
# \  \:\ \:\_\/ /__/:/\:\_\:\ \  \:\    \__\/  \:\/:/ \  \:\ \:\_\/ \__\/  \:\/:/
#  \  \:\_\:\   \__\/  \:\/:/  \  \:\        \__\::/   \  \:\_\:\        \__\::/ 
#   \  \:\/:/        \  \::/    \  \:\       /  /:/     \  \:\/:/        /  /:/  
#    \  \::/          \__\/      \  \:\     /__/:/       \  \::/        /__/:/   
#     \__\/                       \__\/     \__\/         \__\/         \__\/    
# Update 15/04/2020
#=================================================================================================================================

# Please go to here : https://support.microsoft.com/en-gb/help/4026757/windows-10-find-out-how-many-cores-your-processor-has
# to check how many cores your cpu has before you run
number_of_processes=4
#

# For layer 1 2 3 4 5 6 7 8 9 10
similr_max=[1,0.8,0.7,0.6,0.6,0.6,0.6,0.6,0.6,0.6] #<- make sure you input 10 values
C_diff_min=[0,0,0,0,0,0,0,0,0,0]
C_diff_max=[7,7,7,7,7,7,7,7,7,7]
R_diff_min=[0,0,0,0,0,0,0,0,0,0]#<set 1, SPLASH will only choose ring forming template
R_diff_max=[1,1,1,1,1,1,1,1,1,1]
expand_max=[10,10,10,10,10,10,10]

# put SMILES here
input_target=['C=CCN1C2C=C(C(NC3=C4C=CC=C3)C24CC1)C=O']

# Make sure there is not folder with the same name in where you put this script.
folder_name='Biotin'

















#=================================================================================================================================
from multiprocessing import Lock, Process, Queue, current_process, Manager
import time
import queue # imported for using queue.Empty exception
import sys
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
import random
#from keras.models import model_from_json
from rdkit.Chem.Draw import rdMolDraw2D
import random
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import DrawingOptions
import os
import csv
import shutil
import numpy as np
from random import shuffle
import time

def count_ring(smiles):
	try:
		m=Chem.MolFromSmiles(smiles)
		#print('debug smiles',smiles)
		##print('debug m',m)
		ring_count=m.GetRingInfo().NumRings()
		return(ring_count)
	except:
		return(999)
		#-> to make the condition about ring fail,

def count_C(smiles):
	p=smiles
	pc=p.count('C')+p.count('c')-p.count('Cl')-p.count('Cs')-p.count('Cd')-p.count('Cu')-p.count('Co')
	return pc

def compare(smiles_A,smiles_B):
	try:
		mol1 = Chem.MolFromSmiles(smiles_A)
		mol2 = Chem.MolFromSmiles(smiles_B)
		fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1,1)
		fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2,1)
		return DataStructs.TanimotoSimilarity(fp1,fp2)
	except:
		return 0
		#Smiles wrong -> remove this result
def Draw_mol(smiles,filename):
	try:
		m1 = AllChem.MolFromSmiles(smiles)
		mc = rdMolDraw2D.PrepareMolForDrawing(m1)
		drawer = Draw.MolDraw2DCairo(300, 300)
		drawer.DrawMolecule(mc)
		drawer.FinishDrawing()
		output = drawer.GetDrawingText()
		with open(filename,'wb') as pngf:
			pngf.write(output)
	except:
		print('cannot draw from smiles ',smiles)
def Draw_temp(smarts,filename):
	try:
		rxn = AllChem.ReactionFromSmarts(smarts)
		rimage = Draw.ReactionToImage(rxn)
		rimage.save(filename)
	except:
		print('cannot draw from smarts',smarts)
		pass


def apply(temp,target,condition_list,layer):
	similr_max=condition_list[0][layer]
	C_diff_min=condition_list[1][layer] 
	C_diff_max=condition_list[2][layer] 
	R_diff_min=condition_list[3][layer] 
	R_diff_max=condition_list[4][layer]
	#print('condition parameter',similr_max,C_diff_min,C_diff_max,R_diff_min,R_diff_max)
	smarts_input=temp
	nBits=2048
	C_diff_max=7

	rxn = AllChem.ReactionFromSmarts(smarts_input)
	ps = rxn.RunReactants((Chem.MolFromSmiles(target),))
	target_NC=count_C(target)
	target_NR=count_ring(target)
	if ps:

		if len(ps[0])==3:
			p1=(Chem.MolToSmiles(ps[0][0]))
			p2=(Chem.MolToSmiles(ps[0][1]))
			p3=(Chem.MolToSmiles(ps[0][2]))
			output_NC=count_C(p1)+count_C(p2)+count_C(p3)
			output_NR=count_ring(p1)+count_ring(p2)+count_ring(p3)
			if compare(target,p1+'.'+p2+'.'+p3) > similr_max:
				return False
			if not(C_diff_min<np.abs(target_NC - output_NC)<C_diff_max and R_diff_min<(target_NR - output_NR)<R_diff_max and output_NC <= target_NC):
				return False
			else:
				return 3,p1,p2,p3
		if len(ps[0])==2:
			p1=(Chem.MolToSmiles(ps[0][0]))
			p2=(Chem.MolToSmiles(ps[0][1]))
			#print('debug p',p1,p2)
			output_NC=count_C(p1)+count_C(p2)
			output_NR=count_ring(p1)+count_ring(p2)

			#print('from template',smarts,'get reactants',p1,'and',p2)
			if compare(target,p1+'.'+p2) > similr_max:
				return False
			if not(C_diff_min<np.abs(target_NC - output_NC) <C_diff_max and R_diff_min<(target_NR - output_NR)<R_diff_max and output_NC <= target_NC):
				return False
				#print('target',target)
				#print('ouput',p1,p2)
				#print('debug case 2 condition C',C_diff_min<np.abs(target_NC - output_NC) <C_diff_max)
				#print('debug case 2 condition R',R_diff_min<(target_NR - output_NR)<R_diff_max)
				#print('debug case 2 condition C',output_NC <= target_NC)
			else:	
				return 2,p1,p2
		if len(ps[0])==1:
			p1=(Chem.MolToSmiles(ps[0][0]))
			output_NC=count_C(p1)
			output_NR=count_ring(p1)
			#print('from template',smarts,'get reactants',p1)
			if compare(target,p1) > similr_max:
				return False
			if not(C_diff_min<=np.abs(target_NC - output_NC) <C_diff_max and R_diff_min<=(target_NR - output_NR)<R_diff_max and output_NC <= target_NC):
				#print('target',target)
				#print('ouput',p1)
				#print('R_diff_min',R_diff_min,'R_diff_max',R_diff_max)
				#print('target_NR',target_NR,'ourput_NR',output_NR)
				#print('np.abs(target_NR - output_NR)',np.abs(target_NR - output_NR))
				#print('debug case 1 condition C',C_diff_min<np.abs(target_NC - output_NC) <C_diff_max)
				#print('debug case 1 condition R',R_diff_min<(target_NR - output_NR)<R_diff_max)
				#print('debug case 1 condition C',output_NC <= target_NC)
				return False
			else:
				return 1,p1
		else:
			#print('from template',smarts,'cannot get any reactants')
			return False


def plot_n_make_folder(previous_reactant,
						reactants,
						Temp_ID,
						Mol_ID,
						path,
						smarts,
						current_gen_path_list,
						predicted_reaction_list,
						outputimage=True):
	##print('debug 6 path',path)
	##print('debug outputimage assign:',outputimage)
	path_mol_list=[]
	string_cut=20
	if reactants[0]==1:
		##print('debug reactant[0]=1 is True')
		path_mol_1=path+'/'+str(1000000000+Mol_ID)+'_temp_'+str(1000000000+Temp_ID)+'.png'
		path_mol_list.append(path_mol_1[:-string_cut])
		path_temp=path+'/'+'temp_'+str(1000000000+Temp_ID)+'.png'
		if outputimage:
			#print('debug output image is true')
			Draw_mol(reactants[1],path_mol_1)
			#print('debug path mol 1',path_mol_1)
			Draw_temp(smarts[Temp_ID],path_temp)
			predicted_reaction_list.append('Tree:'+'\t'+path+'-'+str(Mol_ID)+'\t'+'temp_ID'+'\t'+str(Temp_ID)+'\t'+previous_reactant+'>>'+reactants[1])
	if reactants[0]==2:
		path_mol_1=path+'/'+str(1000000000+Mol_ID)+'_a_temp_'+str(1000000000+Temp_ID)+'.png'
		path_mol_list.append(path_mol_1[:-string_cut])
		path_mol_2=path+'/'+str(1000000000+Mol_ID)+'_b_temp_'+str(1000000000+Temp_ID)+'.png'
		path_mol_list.append(path_mol_2[:-string_cut])
		path_temp=path+'/'+'temp_'+str(1000000000+Temp_ID)+'.png'
		if outputimage:
			Draw_mol(reactants[1],path_mol_1)
			Draw_mol(reactants[2] ,path_mol_2)
			Draw_temp(smarts[Temp_ID],path_temp)
			predicted_reaction_list.append('Tree:'+'\t'+path+'-'+str(Mol_ID)+'\t'+'temp_ID'+'\t'+str(Temp_ID)+'\t'+previous_reactant+'>>'+reactants[1]+'.'+reactants[2])
	if reactants[0]==3:
		path_mol_1=path+'/'+str(1000000000+Mol_ID)+'_a_temp_'+str(1000000000+Temp_ID)+'.png'
		path_mol_list.append(path_mol_1[:-string_cut])
		path_mol_2=path+'/'+str(1000000000+Mol_ID)+'_b_temp_'+str(1000000000+Temp_ID)+'.png'
		path_mol_list.append(path_mol_2[:-string_cut])
		path_mol_3=path+'/'+str(1000000000+Mol_ID)+'_c_temp_'+str(1000000000+Temp_ID)+'.png'
		path_mol_list.append(path_mol_3[:-string_cut])
		path_temp=path+'/'+'temp_'+str(1000000000+Temp_ID)+'.png'
		if outputimage:
			Draw_mol(reactants[1],path_mol_1)
			Draw_mol(reactants[2] ,path_mol_2)
			Draw_mol(reactants[3] ,path_mol_3)
			Draw_temp(smarts[Temp_ID],path_temp)
			predicted_reaction_list.append('Tree:'+'\t'+path+'-'+str(Mol_ID)+'\t'+'temp_ID'+'\t'+str(Temp_ID)+'\t'+previous_reactant+'>>'+reactants[1]+'.'+reactants[2]+'.'+reactants[3])
	for item in path_mol_list:
		current_gen_path_list.append(item)
def node_expansion(tasks_to_accomplish,
					inputs_q,
					previous_gen_reactants,
					previous_path_list,
					current_gen_reactants,
					current_gen_path_list,
					predicted_reaction_list,
					smarts,
					number_of_processes,
					condition_list,
					layer,
					first_layer):
	loop_limit=condition_list[5][layer]
	#print('start node expansion')
	while True:
		try:
			task = tasks_to_accomplish.get_nowait() #<- task_to_accomplish is queue
			ID=inputs_q.get()
			if len(previous_gen_reactants) > number_of_processes:
				new_ID=ID+number_of_processes
				inputs_q.put(new_ID)
				tasks_to_accomplish.put("Layer "+str(layer+1)+": Task with queque ID " + str(new_ID) +"over total "+str(len(previous_gen_reactants)))
			target=previous_gen_reactants[ID]
			##print('debug 3 previous_path_list',previous_path_list)
			path=previous_path_list[ID]
			#print(target)
			 # <- for later when we need ID for each predicted reactants:
			Draw_mol(target,path+'/000000000.png')
			branch_reactants_no=0
			loop_limit=condition_list[5][layer]
			len_smarts=len(smarts) #<-get number of template
			for Temp_ID in range(0,len_smarts): #<-need template ID for the file name
				if Temp_ID%10000==0 and first_layer:
					os.system('cls' if os.name == 'nt' else 'clear')
					print('first layer will take a while, please wait, especially if the condition you put is too strict or unreasonable')
					print('apply temp ID',Temp_ID,'/',len_smarts)
				temp=smarts[Temp_ID]
				p,r=temp.split('>>')
				if '.' in p: continue
				if branch_reactants_no>=loop_limit:
					break
				out_list=apply(temp,target,condition_list,layer)
				if out_list: #<- check if this molecule appear on the current list
					#print('debug outlist',out_list)
					for item in out_list:
						if item in current_gen_reactants: continue 
					current_gen_reactants.extend(out_list[1:])
					branch_reactants_no+=1
					##print('debug 2 path',path)
					Mol_ID=random.randint(0,999999999)
					plot_n_make_folder(target,
										out_list,
										Temp_ID,
										Mol_ID,
										path,
										smarts,
										current_gen_path_list,
										predicted_reaction_list,
										)

		except queue.Empty:
			break
		else:
			print(task,'are done') # In splash, will print which task for which 
			time.sleep(1)
	return True


def layer_expansion(previous_gen_reactants,previous_path_list,condition_list,number_of_processes,layer,first_layer=False):
	loop_limit=condition_list[5][layer]
	if loop_limit==0:
		sys.exit("No more layer")
	#print('debug 4 previous_path_list:',previous_path_list)
	with open("SMARTS_template_no_chirality.txt", "r") as f:
		smarts=f.readlines()
	tasks_to_accomplish = Queue() #<- this to show which current are running
	inputs_q= Queue()
	processes = []
	manager=Manager()
	current_gen_reactants = manager.list() #<- Make this so that the list can be shared
	current_gen_path_list = manager.list() #<- Make this so that the list can be shared
	predicted_reaction_list = manager.list() #<- Make this so that the list can be shared

	# Append ID for process
	if len(previous_gen_reactants) < number_of_processes:
		for i in range(0,len(previous_gen_reactants)):
			#print('debug 5 previous_gen_reactants',previous_gen_reactants)
			tasks_to_accomplish.put("Layer "+str(layer+1)+": Task with queque ID " + str(i))
			inputs_q.put(i)
	else:
		for i in range(0,number_of_processes):
			# If the lenght of queue is too big, the script will get error
			##print('debug 5 previous_gen_reactants',previous_gen_reactants)
			tasks_to_accomplish.put("Layer "+str(layer+1)+": Task with queque ID " + str(i) +"over total "+str(len(previous_gen_reactants)))
			inputs_q.put(i)
		# creating processes

	for w in range(number_of_processes):
		p = Process(target=node_expansion, args=(tasks_to_accomplish,
													inputs_q,
													previous_gen_reactants,
													previous_path_list,
													current_gen_reactants,
													current_gen_path_list,
													predicted_reaction_list,
													smarts,
													number_of_processes,
													condition_list,
													layer,
													first_layer))
		processes.append(p)
		p.start()
	# completing process
	for p in processes:
		p.join()
	# print the output
	print('current_gen_reactants', current_gen_reactants)
	#print('current_path_list', current_gen_path_list)
	# make folder for next expansion
	for item in current_gen_path_list:
		os.mkdir(item)
	# write reaction path list
	with open(folder_name+'/0_reaction_for_layer'+str(layer+1)+'.txt', 'w') as f:
		for item in predicted_reaction_list:
			f.write("%s\n" % item)
	print('-----------------------------------------------------------------')
	print('-----------------------------END LAYER '+str(layer+1)+'-------------------------')
	print('-----------------------------------------------------------------')
	return  current_gen_reactants,current_gen_path_list,layer+1
if __name__ == '__main__':

	# Get those similarity limit:
	condition_list=[]
	condition_list.extend([similr_max,C_diff_min,C_diff_max,R_diff_min,R_diff_max,expand_max])



	layer=0
	try:
		os.mkdir(folder_name)
	except:
		sys.exit('a Folder with the same name already exist')
	reactants_list,path_list,layer=layer_expansion(input_target,[folder_name],condition_list,number_of_processes,layer,first_layer=True)
	reactants_list,path_list,layer=layer_expansion(reactants_list,path_list,condition_list,number_of_processes,layer)
	reactants_list,path_list,layer=layer_expansion(reactants_list,path_list,condition_list,number_of_processes,layer)
	reactants_list,path_list,layer=layer_expansion(reactants_list,path_list,condition_list,number_of_processes,layer)
	reactants_list,path_list,layer=layer_expansion(reactants_list,path_list,condition_list,number_of_processes,layer)
	reactants_list,path_list,layer=layer_expansion(reactants_list,path_list,condition_list,number_of_processes,layer)
	reactants_list,path_list,layer=layer_expansion(reactants_list,path_list,condition_list,number_of_processes,layer)
	reactants_list,path_list,layer=layer_expansion(reactants_list,path_list,condition_list,number_of_processes,layer)
	reactants_list,path_list,layer=layer_expansion(reactants_list,path_list,condition_list,number_of_processes,layer)
	reactants_list,path_list,layer=layer_expansion(reactants_list,path_list,condition_list,number_of_processes,layer)
