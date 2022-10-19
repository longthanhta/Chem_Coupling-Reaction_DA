import xml.etree.ElementTree as ET
import pandas as pd
from rdkit import Chem


def fix_mol_out_SMILES(mol_lines):
	try:
		lines=mol_lines.split('\n')
	except:
	#print('fail because cannot split file')
		return 'fail because cannot split file'
	#print(lines)
	#print(type(lines))
	del lines[0:4]

	lines.insert(0,'')
	lines.insert(1,'  Reaxys file\n')
	lines.insert(2,'  0  0  0  0  0  0  0  0  0  0999 V3000')
	out='\n'.join(lines)
	#print('mol_lines',out)
	with open('temp.mol', 'w') as the_file:
	   the_file.write(out)
	try:
		m=Chem.MolFromMolFile('temp.mol')
		SMILES=Chem.MolToSmiles(m)
		#print('SMILES',SMILES)
		return SMILES
	except:
		#print('fail because cannot get mol from SMILES')
		return 'fail because cannot get mol from SMILES'




file = open(f"set_1.xml",'r',encoding='utf8')
string = file.read()
for tag in ['<sub>','</sub>','<sup>','</sup>','<i>','</i>','<hi>','</hi>',
	'</SUB>','<SUB>','</SUP>','<SUP>']:
	string = string.replace(tag,'')
string = string.encode('ascii','ignore').decode()
root = ET.fromstring(string)


dict_out={
	'RXN_ID':[],
	'RCT':[],
	'PRD':[],
	'RXD':[],
	'CAT':[],
	'RGT':[],
	'SOL':[],
	'TXT':[],
	'YIELD':[],
	'DOI':[],
	'REF':[],
	}

for reaction in root.iter('reaction'):
	# find reaction ID
	print('\nreaction',reaction.attrib)
	rx_id = reaction.find('RX').find('RX.ID').text
	print('rx_id',rx_id)
	# find reactan and product
	rct = []
	pro = []


	if(reaction.find('RY') is not None):
		for entity in reaction.find('RY'):
			if 'RY.RCT' in entity.tag:
				try: rct.append(entity.text)
				except: rct.append('')
			elif 'RY.PRO' in entity.tag:
				try: pro.append(entity.text)
				except: pro.append('')
				# find different set of condition
		print('number of RXD',len(reaction.findall('RXD')))

	for rc_i in range(len(rct)):
		rct[rc_i] = fix_mol_out_SMILES(rct[rc_i])

	for pr_i in range(len(pro)):
		pro[pr_i] = fix_mol_out_SMILES(pro[pr_i])


	rct_str = '.'.join(rct)
	pro_str = '.'.join(pro)

	for rxd_i,rxd in enumerate(reaction.findall('RXD')): # each reaction
		print('article',rxd_i+1)
		text_desc=''
		rxdyd = 'na'
		rxd_rxdtxt_lst=rxd.findall('RXD.TXT')
		print('rxd_rxdtxt_lst',rxd_rxdtxt_lst)
		if len(rxd_rxdtxt_lst)!=0:
			rxd_rxdtxt=rxd_rxdtxt_lst[0]
			text_desc=rxd_rxdtxt.text
			print('text_desc',text_desc)
		if rxd.find('RXD01') is not None:
			rxd_rxd01_lst=rxd.find('RXD01')
			if(rxd_rxd01_lst.find('RXD.NYD') is not None):
				rxdyd=rxd_rxd01_lst.find('RXD.NYD').text
			elif(rxd_rxd01_lst.find('RXD.YD') is not None):
				rxdyd=rxd_rxd01_lst.find('RXD.YD').text
		else: continue
		if rxd.findall('RXDS01') is not None:
			rxd_rxds01_lst=rxd.findall('RXDS01')
		else: continue
		extr_cond=False
		for rxd_rxds01_i,rxd_rxds01 in enumerate(rxd_rxds01_lst): # each condition
			# Check the "high light":
			print('rxd_rxds01_i',rxd_rxds01_i)
			if(rxd_rxds01_i>= len(rxd_rxd01_lst)):
				break
			rxd_rxd01=rxd_rxd01_lst[rxd_rxds01_i]
			cond_rgt_lst=[]
			cond_cat_lst=[]
			cond_sol_lst=[]
			cond_cat_str = ''
			cond_rgt_str = ''
			cond_sol_str = ''
			rxd_rxds01_rxd03_lst = rxd_rxds01.findall('RXD03')
			#print('number of reagent',len(rxd_rxds01_rxd03_lst))
			for rxd_rxds01_rxd03 in rxd_rxds01_rxd03_lst:
				rxd_rxds01_rxd03_rxdrgt_lst = rxd_rxds01_rxd03.findall('RXD.RGT')
				if len(rxd_rxds01_rxd03_rxdrgt_lst)==1:
					rxd_rxds01_rxd03_rxdrgt=rxd_rxds01_rxd03_rxdrgt_lst[0]
					rgt_match=rxd_rxds01_rxd03_rxdrgt.attrib
					rgt_text=rxd_rxds01_rxd03_rxdrgt.text
					cond_rgt_lst.append(rgt_text)
					if rgt_match: # FOUND!==========================================
						extr_cond=True

			rxd_rxds01_rxd04_lst = rxd_rxds01.findall('RXD04')
			#print('number of catalyst',len(rxd_rxds01_rxd04_lst))
			for rxd_rxds01_rxd04 in rxd_rxds01_rxd04_lst:
				rxd_rxds01_rxd04_rxdcat_lst = rxd_rxds01_rxd04.findall('RXD.CAT')
				if len(rxd_rxds01_rxd04_rxdcat_lst)==1:
					rxd_rxds01_rxd04_rxdcat=rxd_rxds01_rxd04_rxdcat_lst[0]
					cat_match=rxd_rxds01_rxd04_rxdcat.attrib
					cat_text=rxd_rxds01_rxd04_rxdcat.text
					cond_cat_lst.append(cat_text)
					if cat_match: # FOUND!==========================================
						extr_cond=True

			rxd_rxds01_rxd05_lst = rxd_rxds01.findall('RXD05')
			#print('number of solvent',len(rxd_rxds01_rxd05_lst))
			for rxd_rxds01_rxd05 in rxd_rxds01_rxd05_lst:
				rxd_rxds01_rxd05_rxdsol_lst = rxd_rxds01_rxd05.findall('RXD.SOL')
				if len(rxd_rxds01_rxd05_rxdsol_lst)==1:
					rxd_rxds01_rxd05_rxdsol=rxd_rxds01_rxd05_rxdsol_lst[0]
					sol_match=rxd_rxds01_rxd05_rxdsol.attrib
					sol_text=rxd_rxds01_rxd05_rxdsol.text
					cond_sol_lst.append(sol_text)
					if sol_match: # FOUND!==========================================
						extr_cond=True
			#print('contain highligh or not?',extr_cond)

			if extr_cond:
			#	dict_out['EXTR_COND'].append('True')
			#else:
			#	dict_out['EXTR_COND'].append('False')
				cond_rgt_str='.'.join(cond_rgt_lst)
				cond_cat_str='.'.join(cond_cat_lst)
				cond_sol_str='.'.join(cond_sol_lst)
				print('article',rxd_i+1)
				print('_ cond',rxd_rxds01_i+1)
				print('_ _ rgt',cond_rgt_str)
				print('_ _ cat',cond_cat_str)
				print('_ _ sol',cond_sol_str)

				rxd_ref_lst = []
				rxd_doi = 'na'
				if(rxd.find('citations').find('citation')) is not None:
					rxd_citations_cit = rxd.find('citations').find('citation').find('CIT')
					if(rxd_citations_cit.find('CIT.AU') is not None):
						rxd_ref_lst.append(rxd_citations_cit.find('CIT.AU').text)
					if(rxd_citations_cit.find('CIT01') is not None):
						rxd_citations_cit_cit01 = rxd_citations_cit.find('CIT01')
						if(rxd_citations_cit_cit01.find('CIT.JTS') is not None):
							rxd_ref_lst.append(rxd_citations_cit_cit01.find('CIT.JTS').text)
						if(rxd_citations_cit_cit01.find('CIT.LA') is not None):
							rxd_ref_lst.append(rxd_citations_cit_cit01.find('CIT.LA').text)
						if(rxd_citations_cit_cit01.find('CIT.PY') is not None):
							rxd_ref_lst.append(rxd_citations_cit_cit01.find('CIT.PY').text)
						if(rxd_citations_cit_cit01.find('CIT.VL') is not None):
							rxd_ref_lst.append(rxd_citations_cit_cit01.find('CIT.VL').text)
						if(rxd_citations_cit_cit01.find('CIT.NB') is not None):
							rxd_ref_lst.append(rxd_citations_cit_cit01.find('CIT.NB').text)
						if(rxd_citations_cit_cit01.find('CIT.PAG') is not None):
							rxd_ref_lst.append(rxd_citations_cit_cit01.find('CIT.PAG').text)
						if(rxd_citations_cit_cit01.find('CIT.DOI') is not None):
							rxd_doi = rxd_citations_cit_cit01.find('CIT.DOI').text

				rxd_ref = ' '.join(rxd_ref_lst)

				dict_out['DOI'].append(rxd_doi)
				dict_out['RXN_ID'].append(rx_id)
				dict_out['RCT'].append(rct_str)
				dict_out['PRD'].append(pro_str)
				dict_out['CAT'].append(cond_cat_str)
				dict_out['RGT'].append(cond_rgt_str)
				dict_out['SOL'].append(cond_sol_str)
				dict_out['YIELD'].append(rxdyd)
				dict_out['REF'].append(rxd_ref)
				dict_out['RXD'].append(rxd_i+1)
				dict_out['TXT'].append(text_desc)
	#break
df = pd.DataFrame(dict_out)

df.to_excel('output.xlsx')
