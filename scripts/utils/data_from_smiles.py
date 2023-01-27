from rdkit import Chem
from rdkit.Chem import AllChem


def data_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol,radius=2,nBits=2048)
    return fingerprint

