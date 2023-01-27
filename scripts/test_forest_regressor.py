import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from joblib import load
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import PredictionErrorDisplay
from utils import data_from_smiles


def test_forest_regressor(model_path, data_path, out_path):
    con = sqlite3.connect(data_path)
    cur = con.cursor()

    query = """
        SELECT *
        FROM assays a
        WHERE data_split = 'test'
    """
    test = pd.read_sql(query, con, index_col='CID')
    X = [data_from_smiles(s) for s in test['SMILES']]
    y = test['pIC50']
    regr = load(model_path)
    y_pred = regr.predict(X)
    print('Score on test data:')
    print(regr.score(X, y))
    rfr_disp = PredictionErrorDisplay(y_true=y, y_pred=y_pred)
    rfr_disp.plot()
    plt.savefig(out_path)


test_forest_regressor(snakemake.input[0], snakemake.input[1], snakemake.output[0])
