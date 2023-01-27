import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from joblib import load
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import RocCurveDisplay
from utils import data_from_smiles


def test_forest_classifier(model_path, data_path, out_path):
    con = sqlite3.connect(data_path)
    cur = con.cursor()

    query = """
        SELECT *
        FROM assays a
        WHERE data_split = 'test'
    """
    test = pd.read_sql(query, con, index_col='CID')
    X = [data_from_smiles(s) for s in test['SMILES']]
    y = test['did_bind']
    clf = load(model_path)
    print('Score on test data:')
    print(clf.score(X, y))
    rfc_disp = RocCurveDisplay.from_estimator(clf, X, y)
    rfc_disp.figure_.suptitle("ROC curve")
    plt.plot([0, 1], [0, 1], "k--")
    plt.savefig(out_path)


test_forest_classifier(snakemake.input[0], snakemake.input[1], snakemake.output[0])
