import pandas as pd

filename = 'extracted_SMILES_reaxys_NiCat_only_cat_col.csv'
df = pd.read_csv(filename, header = 1,
                 names = ['ID','reaction','PN','DOI','Catalysts'])
print(df)
catalysts = df.Catalysts.str.split('\n')
newCata = []
print(catalysts)
for catas in catalysts:
    print(catas)
    noNickel = []
    for oneCata in catas:
        if (oneCata.find('Ni') == -1 and oneCata.find('nickel') == -1):
            noNickel.append(oneCata)
    for cata in noNickel:
        catas.remove(cata)
    catas = '\n'.join(catas)
    newCata.append(catas)

df = pd.DataFrame({'reaction':df.reaction,
                   'ID':df.ID,
                   'PN':df.PN,
                   'DOI':df.DOI,
                   'catalyst':newCata})
df.to_csv('new_'+filename)
