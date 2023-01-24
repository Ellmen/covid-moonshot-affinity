import sqlite3
import pandas as pd


def explore_data(data_path, out_path):
    con = sqlite3.connect(data_path)
    cur = con.cursor()
    txt = ''

    # All compounds
    query = """
        SELECT count(*)
        FROM compounds
    """
    res = cur.execute(query)
    rows = res.fetchone()
    num_compounds = rows[0]
    txt += f'There are {num_compounds} compounds\n'

    # All compounds for Lipinski's rule of 5
    query = """
        SELECT count(*)
        FROM compounds
        WHERE HBD <= 5 AND HBA <= 10 AND MW <= 500 AND logP <= 5
    """
    res = cur.execute(query)
    rows = res.fetchone()
    num_compounds = rows[0]
    txt += f'There are {num_compounds} compounds which satisfy Lipinksi\'s rule of 5\n'

    # All compounds which violate Lipinski's rule of 5 at most once
    query = """
        SELECT count(*)
        FROM compounds
        WHERE
            CASE WHEN HBD > 5 THEN 1 ELSE 0 END +
            CASE WHEN HBA > 10 THEN 1 ELSE 0 END +
            CASE WHEN MW > 500 THEN 1 ELSE 0 END +
            CASE WHEN logP > 5 THEN 1 ELSE 0 END
        <= 1;
    """
    res = cur.execute(query)
    rows = res.fetchone()
    num_compounds = rows[0]
    txt += f'There are {num_compounds} compounds which violate Lipinksi\'s rule of 5 at most once\n'
    con.close()
    with open(out_path, 'w') as f:
        f.write(txt)

explore_data(snakemake.input[0], snakemake.output[0])
