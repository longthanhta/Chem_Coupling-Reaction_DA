from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

def out(smiles):
    mol = Chem.MolFromSmiles(smiles)
    
    generator = rdFingerprintGenerator.GetMorganGenerator(radius = 1)
    fingerprint = generator.GetSparseCountFingerprint(mol)
    non_zero = fingerprint.GetNonzeroElements()
    
    print(non_zero)


out('CC')
out('CO')
out('CCC')
out('CCO')
out('CO')
out('O=C(O)CCCCC1C(NC(N2)=O)C2CS1')
