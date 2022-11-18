from rdkit import Chem
from rdkit.Chem import AllChem
import random, math, time
import numpy as np

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

### Preprocessing, prepare data, functions
# Get result from RDKI
def test_template(ind,target,neighbor,template): 
	p, r = template[ind].split('>>')
	if '.' in p: # Dont apply template with two product
		return False
	rxn = AllChem.ReactionFromSmarts(template[ind])
	ps = rxn.RunReactants((Chem.MolFromSmiles(target),))
	if ps:
		if len(ps[0])==2:
			p1=(Chem.MolToSmiles(ps[0][0]))
			p2=(Chem.MolToSmiles(ps[0][1]))
            
			# return p1,p2
			# Now searching in list of reactions we collected
			if Get_FP(p1) in neighbor and Get_FP(p2) in neighbor: 
            # TODO, have to check with finger print! DONE, to enhance accuracy
            # smile_to_id[p1] in lib_react[target] and smile_to_id[p2] in lib_react[target]:
				return True
			else: 
				return False

		if len(ps[0])==1:
			p1=(Chem.MolToSmiles(ps[0][0]))
			# return p1
			if Get_FP(p1) in neighbor:
				return True
			else:
				return False
		else:
			return False
        

def Get_CT(smile,chars,id_to_efp,smile_to_id,lib_react,template): # Check Template (CT)
    neighbor = [id_to_efp[id_] for id_ in lib_react[smile_to_id[smile]]]
    score = []
    for j in range(NUM_OF_TEMPLATE):
        try:
            outcome = test_template(j,smile,neighbor,template)
        except:
            outcome = False
        if outcome:
            score.append(1)
        else:
            score.append(0)
    return encrypt_ct(score,chars)
    
    
def Get_FP(smile,chars,nBits=2048):
    return encrypt_fp(np.array(AllChem.GetMorganFingerprintAsBitVect(            AllChem.MolFromSmiles(smile), 2, nBits,useChirality=False), dtype=np.byte),chars)


### Encryption
def six_bits_to_code(x:iter)->int:
    # size of x must be 6, 2^6 = 64 characters
    # note that the order of conversion is quite reversed
    return sum([x[i]*pow(2,i) for i in range(6)])


def encrypt_fp(fp:iter,chars:str)->str:
    
    """ Convert binary string to general string {0,1} -> {a:z,A:Z,0:9,$@} 
    (64 characters). Since we have 2048 bits, we ADD 4 0's at the end
    to take 342 6-bit phrases. When use CONVERT BACK to fp, PLEASE REMEMBER 
    to REMOVE 4 0's at the end"""
    
    # fp must have size 2048
    fp = list(fp)
    fp += [0,0,0,0]
    out = ""
    for i in range(342):
        out += chars[six_bits_to_code(fp[i*6:i*6+6])]
    return out


def encrypt_ct(score:iter,chars:str)->str:
    
    # score must have size 193538 = NUM_OF_TEMPLATE
    score = list(score)
    score += [0,0,0,0]
    out = ""
    for i in range(32257):
        out += chars[six_bits_to_code(score[i*6:i*6+6])]
    return out


### Task
def task(list_smile_target,chars,id_to_efp,smile_to_id,
         lib_react,template,list_output):
    for smile_ in list_smile_target:
        list_output.append([smile_to_id[smile_],id_to_efp[smile_to_id[smile_]],Get_CT(smile_,chars,id_to_efp,smile_to_id,lib_react,template)])
        

if __name__ == "__main__":
    
    t = time.time()
    # Import template data
    with open("SMARTS_template.txt", "r") as f: # TODO
    	template=f.readlines() # so we don't need to remove '\n'
    NUM_OF_TEMPLATE = len(template)
    
    # Import lib!
    # Dictionay from chem_id to SMARTS and vice versa
    dict_file = open("new_chem_dictionary.txt",'r')
    smile_to_id = {}
    id_to_smile = {}
    NUM_CHEM = 0
    for line in dict_file: 
        NUM_CHEM += 1
        [id_ ,smile_] = line[:-1].split(',') 	# eliminate '\n', split
        smile_to_id[smile_] = id_   # from smile to id_
        id_to_smile[id_] = smile_   # from id_ to smile
    dict_file.close()
    
    # Reaction in simple form
    reaction_file = open("new_reaction_reverse.txt",'r')
    lib_react = {}
    for line in reaction_file:
        line = line[:-1]        # eliminate '\n'
        [pro, rea] = line.split(':')
        lib_react[pro] = rea.split(',')
        while '' in lib_react[pro]:
            lib_react[pro].remove('')
    reaction_file.close()
    
    # Finger print
    finger_print_file = open("encrypted_id_fp.csv",'r')
    id_to_efp = dict()
    for line in finger_print_file:
        [id_ ,efp_] = line[:-1].split(',')
        id_to_efp[id_] = efp_
    finger_print_file.close()

    # Take training target and checking target
    TRAIN_SIZE = 2
    TRAIN_CLUSTER_SIZE = math.floor(NUM_CHEM/TRAIN_SIZE)
    train_target = [id_to_smile[str(random.randint(0,TRAIN_CLUSTER_SIZE-1)+i*TRAIN_CLUSTER_SIZE)] for i in range(TRAIN_SIZE)]
    random.shuffle(train_target)
    
    # Preparation
    chars = "abcdefghijklnmopqrstuvwxyzABCDEFJHIJKLNMOPQRSTUVWXYZ0123456789$@"
    id_X_Y = list()
    task(train_target,chars,id_to_efp,smile_to_id,lib_react,template,id_X_Y)
    out = open("data_id_X_Y.csv",'w')
    for triple in id_X_Y:
        out.write(','.join(triple)+'\n')
    out.close()
    print("Finish in:",time.time()-t)

""" I write code in the setting for multiprocessing, say all tasks can be seperated from each other and the main process!


"""

