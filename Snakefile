rule all:
    input:
        "outputs/lipinski.txt"

rule build_db:
    input:
        "data/raw/activity_data.csv"
    output:
        "data/sqlite/activity_data_raw.db"
    script:
        "scripts/build_db.py"

rule add_descriptors:
    input:
        "data/sqlite/activity_data_raw.db"
    output:
        "data/sqlite/activity_data.db"
    script:
        "scripts/add_descriptors.py"

rule explore_data:
    input:
        "data/sqlite/activity_data.db"
    output:
        "outputs/lipinski.txt"
    script:
        "scripts/explore_data.py"

rule plot_regression:
    input:
        "models/linear_model.joblib",
        "data/cleaned/airquality_cleaned.csv"
    output:
        "outputs/reg_plot.pdf"
    script:
        "scripts/plot_linear_model.py"
