# this folder is for extracting Ar6nonHet-Ar6nonHet reaction data and add information
raw_data 
->
1. filter_HasAr6_AND_add_SubInfo_HetStr.py
        - filter raw_data to only contain reaction with 6-member aromatic ring in either reactant
        - then add subsitutent information 
        - and add a string description of Heteroatom in the aromatic ring
-> 
2. filter_Ar6nonHet-Ar6nonHet.py
        - filter data to only contain reaction with Phenyl-Phenyl coupling
-> 
3. add_Smil2SubFeat.py
        - add substituent features from smiles
