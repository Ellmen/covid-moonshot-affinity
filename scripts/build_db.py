import sqlite3
import pandas as pd


def build_db(data_path, out_path):
    con = sqlite3.connect(out_path)
    cur = con.cursor()
    df = pd.read_csv(data_path, index_col='CID')
    assays_df = df.drop(columns=['SMILES'])
    compounds_df = df[['SMILES']]
    assays_table = 'assays'
    assays_df.to_sql(assays_table, con, if_exists='replace', index=True)
    compounds_table = 'compounds'
    compounds_df.to_sql(compounds_table, con, if_exists='replace', index=True)
    con.close()


build_db(snakemake.input[0], snakemake.output[0])
