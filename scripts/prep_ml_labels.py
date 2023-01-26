import shutil

import sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors


def add_descriptors(data_path, out_path):
    shutil.copyfile(data_path, out_path)
    con = sqlite3.connect(out_path)
    cur = con.cursor()
    assays_df = pd.read_sql('SELECT * FROM assays', con, index_col='CID')
    assays_df['data_split'] = None
    assays_df['did_bind'] = 0
    assays_df = assays_df[assays_df['f_avg_IC50'].notnull()]

    for index, row in assays_df.iterrows():
        assays_df.at[index, 'pIC50'] = -np.log10(row['f_avg_IC50'])
        assays_df.at[index, 'did_bind'] = 0 if row['f_avg_IC50'] >= 99 else 1

    train, validate, test = np.split(
        assays_df.sample(frac=1, random_state=42),
        [int(0.8*len(assays_df)), int(0.9*len(assays_df))]
    )
    train['data_split'] = 'train'
    validate['data_split'] = 'validate'
    test['data_split'] = 'test'
    assays_df = pd.concat([train, validate, test])
    assays_df.to_sql('assays', con, if_exists='replace', index=True)
    con.close()


add_descriptors(snakemake.input[0], snakemake.output[0])
