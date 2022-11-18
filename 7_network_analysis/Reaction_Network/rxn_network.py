from rdkit.Chem import AllChem
import networkx as nx
import pandas as pd

def TanimotoSimilarity(fp1, fp2):
	#print('fp1',fp1)
	#print('fp2',fp2)
	if(len(fp1) != len(fp2)):
		print("Finger prints have different length")
		return 0
	intersection = [fp1[i]*fp2[i] for i in range(len(fp1))]
	sum_intersection = sum(intersection)
	sum_union = sum(fp1)+sum(fp2)-sum_intersection
	if sum_union == 0: 
		# that means all entries of fp1 and fp2 are 0
		return 1
	else:
		tanimoto = sum_intersection/sum_union
		return tanimoto

def MorganFingerPrint(smile,radius=2,nBits=2048):
	mol =AllChem.MolFromSmiles(smile)
	mol.UpdatePropertyCache()
	fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits = nBits)
	return fp

rxn_list = pd.read_excel('7_Aug_labeled_N_to_Ni_bond_extraction.xlsx').fillna('')

edge_list = [] # paired up by reaxys, directed
similarity_threshold = 0.8
in_degree = []
out_degree = []

for rxn1 in rxn_list.index:
	for rxn2 in rxn_list.index:
		# check whether a product of rxn 1 is a similar to rct rxn 2
		# product of rxn 1
		print('Reaction 1:',rxn_list['SMILES'][rxn1])
		prd1 = rxn_list['SMILES'][rxn1].split('>>')[1].split('.')
		prd1 = [MorganFingerPrint(s) for s in prd1]

		# reactants of rxn 2
		print('Reaction 2:',rxn_list['SMILES'][rxn2])
		rct2 = rxn_list['SMILES'][rxn2].split('>>')[0].split('.')
		rct2 = [MorganFingerPrint(s) for s in rct2]

		# check similarity
		for p in prd1:
			for r in rct2:
				if TanimotoSimilarity(p,r) > similarity_threshold:
					# if rxn1 has a product similar to a reactant of rxn2
					# we lay a directed edge from rxn1 to rxn2, with node label as Reaxys ID
					edge_list.append((rxn_list['RXD'][rxn1],rxn_list['RXD'][rxn2]))
					break

# save edge list, so do not need to run the above part again
out = open('edge_list.csv','w')
for edge in edge_list:
	out.write(f'{edge[0]},{edge[1]}\n')
out.close()

# build graph
DG = nx.DiGraph()
DG.add_edge_from(edge_list)
options = {
    'node_color': 'blue',
    'node_size': 1200,
    'width': 1,
    'with_labels': True
}
nx.draw_spectral(DG, **options)
plt.savefig('reaction_network.png')