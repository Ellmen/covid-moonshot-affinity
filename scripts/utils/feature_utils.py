# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 10:30:52 2023

@author: HP ZBOOK
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.feature_selection import RFE, SelectKBest, f_classif
from sklearn.svm import SVR
from sklearn.linear_model import LogisticRegression

def compute_descriptors(molecule):
    descriptors ={}
    for d in Descriptors.descList:
        try:
            value = d[1](molecule)
            descriptors[d[0]] = value

        except:
            value = np.nan
            descriptors[d[0]] = np.nan
    descriptors = pd.Series(descriptors)
    return descriptors

def pack_df(df_with_SMILES):
    descriptorsdf = pd.DataFrame()
    for index, row in df_with_SMILES.iterrows():
    # Get the SMILES string of the compound
        smiles = row['SMILES']
        # Create a rdkit molecule object from the SMILES string
        mol = Chem.MolFromSmiles(smiles)
        des = compute_descriptors(mol)
        for i in range(len(des)):
                descriptorsdf.at[index,'{}'.format(des.index[i])] = des.values[i]
    return descriptorsdf

def feature_selection(Descriptors_df, y_df, k_feature = 10, task = 'regression', method = 'simple', rfe_step = 0.2):
    if task == 'regression':
        estimator = SVR(kernel="linear")
    if task == 'classification':
        estimator = LogisticRegression()
    if method == 'simple':
        selector = SelectKBest(f_classif, k=k_feature)
    if method == 'rfe':
        selector = RFE(estimator, n_features_to_select=k_feature, step=rfe_step)

    selector.fit(Descriptors_df, y_df)
    # Get the indices of the top k features
    top_k_idx = selector.get_support(indices=True)

    print("Top", k_feature, "features:", Descriptors_df.columns[top_k_idx])

    return Descriptors_df.columns[top_k_idx]


