import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Function to compute pIC50 values
def compute_pIC50(ic50):
    if pd.isnull(ic50):
        return np.nan
    else:
        return -np.log10(ic50)

def plot_ic50(data_path, out_path):
    con = sqlite3.connect(data_path)
    cur = con.cursor()

    # Select IC50 values from the assays table
    query = """
        SELECT f_avg_IC50
        FROM assays
    """
    res = cur.execute(query)
    ic50_values = [row[0] for row in res.fetchall()]

    # Compute pIC50 values
    pic50_values = [compute_pIC50(item) for item in ic50_values]

    # replace inf values with NaN
    pic50_values = np.where(np.isinf(pic50_values), np.nan, pic50_values)

    plt.hist(pic50_values, bins = 20)
    plt.xlabel('pIC50')
    plt.ylabel('Frequency')
    plt.savefig(out_path)


plot_ic50(snakemake.input[0], snakemake.output[0])
