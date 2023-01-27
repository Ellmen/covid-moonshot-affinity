import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from joblib import dump
from sklearn.ensemble import RandomForestRegressor
from utils import data_from_smiles


def forest_regressor(data_path, out_path, model_out_path):
    con = sqlite3.connect(data_path)
    cur = con.cursor()

    query = """
        SELECT *
        FROM assays a
        WHERE data_split = 'train'
    """
    train = pd.read_sql(query, con, index_col='CID')
    query = """
        SELECT *
        FROM assays
        WHERE data_split = 'validate'
    """
    validate = pd.read_sql(query, con, index_col='CID')
    X = [data_from_smiles(s) for s in train['SMILES']]
    y = train['pIC50']
    Xv = [data_from_smiles(s) for s in validate['SMILES']]
    yv = validate['pIC50']
    best_val_score = 0
    best_model = None
    val_scores = []
    max_depths = []
    for max_depth in range(1, 11):
        print(max_depth)
        regr = RandomForestRegressor(max_depth=max_depth, random_state=0)
        regr.fit(X, y)
        print('Score on training data:')
        print(regr.score(X, y))
        print('Score on validation data:')
        val_score = regr.score(Xv, yv)
        print(val_score)
        val_scores.append(val_score)
        max_depths.append(max_depth)
        if val_score > best_val_score:
            best_val_score = val_score
            best_model = regr
    plt.plot(max_depths, val_scores)
    plt.ylabel('Mean validation score')
    plt.xlabel('Max depth')
    plt.savefig(out_path)
    dump(best_model, model_out_path)


forest_regressor(snakemake.input[0], snakemake.output[0], snakemake.output[1])
