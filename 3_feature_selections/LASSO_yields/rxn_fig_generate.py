import pandas as pd
import os
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from PIL import Image
from PIL import Image, ImageDraw, ImageFont

def defaultDrawOptions():
    '''This function returns an RDKit drawing options object with
    default drawing options.'''

    opts = Draw.DrawingOptions()
    # opts.elemDict = defaultdict(lambda: (0,0,0)) # all atoms are black
    opts.noCarbonSymbols = True
    opts.selectColor = (1, 1, 0)
    opts.wedgeBonds = True
    return opts


def draw_reaction(input_smiles,rxn_path):
    print('input_smiles',input_smiles)
    rxn = AllChem.ReactionFromSmarts(input_smiles,useSmiles=True)
    d2d = Draw.MolDraw2DCairo(800,200)
    d2d.drawOptions().useBWAtomPalette()
    options=defaultDrawOptions()
    d2d.DrawReaction(rxn,highlightByReactant=True)
    png = d2d.GetDrawingText()
    open(rxn_path,'wb+').write(png)
    print('ploted ','in',rxn_path)
    #add_text(EG+'  ('+row['AUT']+'  '+row['REF']+')',img_path)

def add_text(full_text,path_img):
    img = Image.open(path_img)
    d1 = ImageDraw.Draw(img)
    myFont = ImageFont.truetype('Aaargh.ttf', 40)
    d1.text((0, 0), "Leaving group: "+full_text, font=myFont)
    img.save(path_img)

data_df=pd.read_excel('CC_CC_data_Ph_Ph_yield_sub.xlsx')
print(data_df)

try:
    os.mkdir('1_reactions')
except:
    print('folder reacetion already exists')

for index,row in data_df.iterrows():
    AMSMILES=row['AAM']
    rxn_path=str(index)
    draw_reaction(AMSMILES,'1_reactions/'+rxn_path+'.png')
    add_text(row['FILENAME'],'1_reactions/'+rxn_path+'.png')
