import pandas as pd
import os
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole

def defaultDrawOptions():
    '''This function returns an RDKit drawing options object with
    default drawing options.'''

    opts = Draw.DrawingOptions()
    # opts.elemDict = defaultdict(lambda: (0,0,0)) # all atoms are black
    opts.noCarbonSymbols = True
    opts.selectColor = (1, 0, 0)
    opts.wedgeBonds = True
    return opts


def draw_reaction(input_smiles,rxn_path):
    print('input_smiles',input_smiles)
    rxn = AllChem.ReactionFromSmarts(input_smiles,useSmiles=True)
    d2d = Draw.MolDraw2DCairo(1600,300)
    d2d.drawOptions().useBWAtomPalette()
    options=defaultDrawOptions()
    d2d.DrawReaction(rxn,highlightByReactant=True)
    png = d2d.GetDrawingText()
    open(rxn_path,'wb+').write(png)
    print('ploted ','in',rxn_path)
    #add_text(EG+'  ('+row['AUT']+'  '+row['REF']+')',img_path)



data_df=pd.read_excel('all_data_2020_Jan_1.xlsx',sheet_name='Identified CC')
print(data_df)

try:
    os.mkdir('1_reactions')
except:
    print('folder reacetion already exists')

for index,row in data_df.iterrows():
    AMSMILES=row['AAM']
    rxn_path=str(row['OID'])
    draw_reaction(AMSMILES,'1_reactions/'+rxn_path+'.png')
