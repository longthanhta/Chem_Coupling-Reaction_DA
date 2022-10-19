import re
from numpy.random import shuffle
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import ChiralType
from IPython.display import SVG
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw
import numpy as np
import pandas as pd

def Draw_mol(smiles):
    m1 = AllChem.MolFromSmiles(smiles)
    mc = rdMolDraw2D.PrepareMolForDrawing(m1)
    #drawer = Draw.MolDraw2DCairo(300, 300)
    drawer = rdMolDraw2D.MolDraw2DSVG(300,300)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # fix the svg string and display it
    display(SVG(svg.replace('svg:','')))
    

def Get_FP(smiles,nBits=8):
    fp=(np.array(AllChem.GetMorganFingerprintAsBitVect(
    AllChem.MolFromSmiles(smiles), 1
    , nBits,useChirality=False), dtype=np.int))
    fp=np.array(fp)
    Draw_mol(smiles)
    lst=[]
    for i in range(0,len(fp)):
        if fp[i] == 1:
            lst.append(i+1)
    lst.append(smiles)
    return fp

FP1=Get_FP('CC(C)C1=CN=CC=C1')
