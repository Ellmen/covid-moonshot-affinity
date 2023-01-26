import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_ic50(data_path, out_path):
    con = sqlite3.connect(data_path)
    cur = con.cursor()

    # Select IC50 values from the assays table
    query = """
        SELECT pIC50
        FROM assays
        WHERE did_bind = 1
    """
    res = cur.execute(query)
    pic50_values = [row[0] for row in res.fetchall()]

    plt.hist(pic50_values, bins = 20)
    plt.xlabel('pIC50')
    plt.ylabel('Frequency')
    plt.savefig(out_path)


plot_ic50(snakemake.input[0], snakemake.output[0])
