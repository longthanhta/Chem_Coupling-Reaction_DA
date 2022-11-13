from rdkit import Chem
from rdkit.Chem import AllChem
from keras import Sequential
from keras.layers import Dense
import random, math, time
import numpy as np

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

### Preprocessing, prepare data, functions
# Get result from RDKI
def Get(ind,target): 
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
			if smile_to_id[p1] in lib_react[target] and smile_to_id[p2] in lib_react[target]: # TODO
				return True
			else: 
				return False

		if len(ps[0])==1:
			p1=(Chem.MolToSmiles(ps[0][0]))
			# return p1
			if smile_to_id[p1] in lib_react[target]:
				return True
			else:
				return False
		else:
			return False

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

# Take training target and checking target
TRAIN_SIZE = 1000
CHECK_SIZE = 500
TRAIN_CLUSTER_SIZE = math.floor(NUM_CHEM/TRAIN_SIZE)
CHECK_CLUSTER_SIZE = math.floor(NUM_CHEM/CHECK_SIZE)
train_target = [str(random.randint(0,TRAIN_CLUSTER_SIZE-1)+i*TRAIN_CLUSTER_SIZE) for i in range(TRAIN_SIZE)]
check_target = [str(random.randint(0,CHECK_CLUSTER_SIZE-1)+i*CHECK_CLUSTER_SIZE) for i in range(CHECK_SIZE)]
random.shuffle(train_target)
random.shuffle(check_target) # TODO, if check and train overlap??



# Analyze template
#### Generate train X, Y
FINGER_PRINT_SIZE = 2048
# Guarantee that if smile_ is not included in X, so does its score in Y
## Train
fp_lst_train = []
score_train = []
for i in range(TRAIN_SIZE):
	print('====== Generating data with training target no',i,'=======')
	try: 
		fp_lst_train.append(np.array(AllChem.GetMorganFingerprintAsBitVect(
			AllChem.MolFromSmiles(id_to_smile[train_target[i]]), 2, nBits=FINGER_PRINT_SIZE,useChirality=False), dtype=np.byte)) # TODO, load it from encrypted_id_smile
		smile = id_to_smile[train_target[i]]
		score = []
		for j in range(NUM_OF_TEMPLATE):
			if j % 10000 == 0:
				print(j,'/',NUM_OF_TEMPLATE) # What does this serve

			try:
				outcome = Get(j,input)
			except:
				outcome = False
			if outcome:
				score.append(1)
			else:
				score.append(0)
		score_train.append(score)
	except:
		print(id_to_smile[train_target[i]])
		
X_train = np.array(fp_lst_train)
Y_train = np.array(score_train)	
	

## Checking
fp_lst_check=[]
score_check = []
for i in range(CHECK_SIZE):
	print('====== Generating data with checking target no',i,'=======')
	try: 
		fp_lst_check.append(np.array(AllChem.GetMorganFingerprintAsBitVect(
			AllChem.MolFromSmiles(id_to_smile[check_target[i]]), 2, nBits=FINGER_PRINT_SIZE,useChirality=False), dtype=np.byte))
		smile = id_to_smile[check_target[i]]
		score = []
		for j in range(NUM_OF_TEMPLATE):
			if j % 10000 == 0:
				print(j,'/',NUM_OF_TEMPLATE) # What does this serve

			try:
				outcome = Get(j,input)
			except:
				outcome = False
			if outcome:
				score.append(1)
			else:
				score.append(0)
		score_check.append(score)
	except:
		print(id_to_smile[check_target[i]])
		
X_check = np.array(fp_lst_check)
Y_check = np.array(score_check)	
run_no = 1




#### Neural Network starting from here
# Model
classifier = Sequential()
# Input layer + First Hidden Layer
classifier.add(Dense(300, activation='relu', kernel_initializer='random_normal', input_dim=FINGER_PRINT_SIZE))
# Second  Hidden Layer
classifier.add(Dense(300, activation='relu', kernel_initializer='random_normal'))
# Output Layer
classifier.add(Dense(NUM_OF_TEMPLATE, activation='relu', kernel_initializer='random_normal'))
classifier.compile(optimizer ='adam',loss='binary_crossentropy', metrics =['accuracy'])
classifier.fit(X_train,Y_train, batch_size=100, epochs=20)


### Save model, do prediction

model_json = classifier.to_json()
with open("model.json", "w") as json_file:
		json_file.write(model_json)
classifier.save_weights("model.h5")
print("Saved model to disk")
Y_check_pred=classifier.predict(X_check)

# Convert floating result to binary outcome
c_r = 0
threashold = np.mean(Y_check_pred) # TODO
for item in Y_check_pred:
	c_r += 1
	print('Convert to binary, test result number ',c_r)
	for j in range(len(item)):
		if item[j] > threashold:
			item[j]=1
		else:
			item[j]=0


### Error Analysis

pos_true=[]
pos_false=[]
neg_true=[]
neg_false=[]
num_wrong_prediction=[]
all_false_pred=[]

for i in range(CHECK_SIZE):
	print('======= Testing with target no',i,'========')

	count_pos_true=0
	count_pos_false=0
	count_neg_false=0
	count_neg_true=0
	fail_pred=[]

	predict = Y_check_pred[i]
	real = Y_check[i]

	for k in range(NUM_OF_TEMPLATE):
		pred_score = predict[k]
		real_score = real[k]

		if pred_score==1 and real_score==1:
			count_pos_true+=1
		elif pred_score==1 and real_score==0:
			count_pos_false+=1
			fail_pred.append(str(k))
		elif pred_score==0 and real_score==1:
			count_neg_true+=1
			fail_pred.append(str(k))
		elif pred_score==0 and real_score==0:
			count_neg_false+=1

	num_wrong_prediction.append(count_pos_false+count_neg_true)
	pos_true.append(str(count_pos_true))
	pos_false.append(str(count_pos_false))
	neg_true.append(str(count_neg_true))
	neg_false.append(str(count_neg_false))
	all_false_pred.append(fail_pred)


out = open("summary.csv",'w')
out.write(f'Id,SMILE,Number of wrong predictions (over {NUM_OF_TEMPLATE}),List of fail predictions,"Confusion Matrix (PT,PF,NT,NF)"\n')
for i in range(CHECK_SIZE):
	out.write(check_target[i]+','+id_to_smile[check_target[i]]+','+str(num_wrong_prediction[i])+',"'+','.join(all_false_pred[i])+'","'+pos_true[i]+','+pos_false[i]+','+neg_true[i]+','+neg_false[i]+'"\n')
out.close()
print("Finish:",time.time()-t)


"""
* TRAIN_SIZE = 1000, CHECK_SIZE = 500
Non-multi-tasking: runtime = 22556.057943105698s ~ 6.25h
?? Which task can we apply multithreading, multiprocessing to accelerate

"""