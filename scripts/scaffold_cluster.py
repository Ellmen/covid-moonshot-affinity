import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA, NMF
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold


def ClusterFps(fps,cutoff=0.2):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs


def scaffold_cluster(data_path, out_path):
    con = sqlite3.connect(data_path)
    cur = con.cursor()

    compounds_df = pd.read_sql('SELECT * FROM COMPOUNDS', con, index_col='CID')
    compounds_df['fingerprint'] = None
    for index, row in compounds_df.iterrows():
        smiles = row['SMILES']
        mol = Chem.MolFromSmiles(smiles)
        core = MurckoScaffold.GetScaffoldForMol(mol)
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(core,radius=2,nBits=2048)
        compounds_df.at[index, 'fingerprint'] = fingerprint
    fps = compounds_df['fingerprint']
    num_clusters = []
    cutoffs = []
    for i in range(10):
        cutoff = (i+1)/10
        print(f"Clustering with cutoff {cutoff}")
        cs = ClusterFps(fps, cutoff)
        n_cs = len(cs)
        num_clusters.append(n_cs)
        cutoffs.append(cutoff)
        x = [len(c) for c in cs]
        counts, bins = np.histogram(x)
        plt.pie(x)
        plt.title(
            "Size distribution of clusters\n"
            f"{n_cs} clusters generated with a cutoff of {cutoff}"
        )
        plt.savefig(f'outputs/scaffold_clusters_{cutoff}.pdf')
    plt.clf()
    plt.plot(cutoffs, num_clusters)
    plt.ylabel('Number of clusters')
    plt.xlabel('Cutoff value')
    plt.savefig('outputs/scaffold_cluster_cutoffs.pdf')
    con.close()

scaffold_cluster(snakemake.input[0], snakemake.output[0])
