rule all:
    input:
        "outputs/lipinski.txt",
        "outputs/fingerprint_cluster_cutoffs.pdf",
        "outputs/scaffold_cluster_cutoffs.pdf",
        "outputs/dag.svg"

rule render_dag:
    input:
        "Snakefile"
    output:
        "outputs/dag.svg"
    shell:
        "snakemake --dag | dot -Tsvg > outputs/dag.svg"

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

rule fingerprint_cluster:
    input:
        "data/sqlite/activity_data.db"
    output:
        "outputs/fingerprint_cluster_cutoffs.pdf"
    script:
        "scripts/fingerprint_cluster.py"

rule scaffold_cluster:
    input:
        "data/sqlite/activity_data.db"
    output:
        "outputs/scaffold_cluster_cutoffs.pdf"
    script:
        "scripts/scaffold_cluster.py"

rule plot_regression:
    input:
        "models/linear_model.joblib",
        "data/cleaned/airquality_cleaned.csv"
    output:
        "outputs/reg_plot.pdf"
    script:
        "scripts/plot_linear_model.py"
