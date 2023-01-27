from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors


def compute_best_descriptors(molecule):
    li = ['PEOE_VSA8', 'fr_N_O', 'fr_diazo', 'fr_dihydropyridine', 'fr_hdrzine','fr_ketone_Topliss', 'fr_lactam', 'fr_oxazole', 'fr_sulfonamd','fr_sulfone']
    featurelist = []
    for i, d in enumerate(Descriptors.descList):
        if d[0] in li or True:
            value = d[1](molecule)
            featurelist.append(value)
    return featurelist

def data_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fingerprint = list(AllChem.GetMorganFingerprintAsBitVect(mol,radius=2,nBits=2048))
    top_features = compute_best_descriptors(mol)
    return fingerprint + top_features

