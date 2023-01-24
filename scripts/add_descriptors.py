import shutil

import sqlite3
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors


def add_descriptors(data_path, out_path):
    shutil.copyfile(data_path, out_path)
    con = sqlite3.connect(out_path)
    cur = con.cursor()
    compounds_df = pd.read_sql('SELECT * FROM COMPOUNDS', con, index_col='CID')
    # Create new columns for the molecular descriptors
    compounds_df['HBD'] = None
    compounds_df['HBA'] = None
    compounds_df['MW'] = None
    compounds_df['logP'] = None

    # Iterate over the compounds and calculate the molecular descriptors
    for index, row in compounds_df.iterrows():
        # Get the SMILES string of the compound
        smiles = row['SMILES']
        # Create a rdkit molecule object from the SMILES string
        mol = Chem.MolFromSmiles(smiles)
        # Calculate the number of hydrogen bond donors
        h_bond_donors = Chem.Lipinski.NumHDonors(mol)
        # Calculate the number of hydrogen bond acceptors
        h_bond_acceptors = Chem.Lipinski.NumHAcceptors(mol)
        # Calculate the molecular weight
        mw = Descriptors.MolWt(mol)
        # Calculate the octanol-water partition coefficient
        logp = Descriptors.MolLogP(mol)
        # Add the molecular descriptors to the corresponding columns in the dataframe
        compounds_df.at[index, 'HBD'] = h_bond_donors
        compounds_df.at[index, 'HBA'] = h_bond_acceptors
        compounds_df.at[index, 'MW'] = mw
        compounds_df.at[index, 'logP'] = logp

    compounds_df.to_sql('compounds', con, if_exists='replace', index=True)
    con.close()


add_descriptors(snakemake.input[0], snakemake.output[0])
