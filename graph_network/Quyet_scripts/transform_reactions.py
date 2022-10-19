import pandas as pd

filename = "extracted_SMILES_reaxys.csv"
df = pd.read_csv(filename,dtype = str)
rea_list = []
pro_list = []

for row in df[["ID","reaction","Catalysts"]].iterrows():
    rxn_list.append([row[1]['ID'], row[1]['reaction']
                    ,row[1]['Catalysts']])
## get a list of reaction with id and catalyst imformation

rxn = []
SMILES = []
catalysts = []
final = pd.DataFrame({"reaction":rxn, "catalysts":catalysts, "SMILES":SMILES})
final.to_csv("transformed_reactions.csv")
## create an empty file

cnt = 0
for [rxn_id1, rxn1, cata1] in rxn_list:
    for [rxn_id2, rxn2 , cata2] in rxn_list:
        if (cata2.find('Ni') != -1 or cata2.find('nickel') != -1 or
            cata1.find('Ni') != -1 or cata1.find('nickel') != -1):
            rea = rxn1.split(">>")[0]
            pro = rxn2.split(">>")[1]
            if rea == pro: ## if the reactant is also a product
                rxn.append(rxn_id2 + ">>" + rxn_id1)
                catalysts.append(cata2 + ">>" + cata1)
                SMILES.append(rxn2 + ">>" + rea)
                cnt += 1
                if (cnt % 100 == 0):
                    final = pd.DataFrame({"reaction":rxn, 
                                          "catalysts":catalysts
                                          "SMILES":SMILES})
                    final.to_csv("transformed_reactions.csv", 
                                 mode = 'a',header = False)
                    rxn = []
                    catalysts = []
                    SMILES = []
                    ## write to file for every 100 lines
            
            
final = pd.DataFrame({"reaction":rxn, "catalysts":catalysts,"SMILES":SMILES})
final.to_csv("transformed_reactions.csv", mode = 'a',header = False)
## write to file for the rest lines