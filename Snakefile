rule all:
    input:
        "outputs/lipinski.txt",
        "outputs/fingerprint_cluster_cutoffs.pdf",
        "outputs/scaffold_cluster_cutoffs.pdf",
        "outputs/ic50_distribution.pdf",
        "outputs/ic50_binders_distribution.pdf",
        "outputs/dag.svg",
        "outputs/test_classification_plot.pdf",
        "outputs/test_regression_plot.pdf"

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

rule plot_ic50:
    input:
        "data/sqlite/activity_data_labelled.db"
    output:
        "outputs/ic50_distribution.pdf"
    script:
        "scripts/plot_ic50.py"

rule plot_ic50_binders:
    input:
        "data/sqlite/activity_data_labelled.db"
    output:
        "outputs/ic50_binders_distribution.pdf"
    script:
        "scripts/plot_ic50_binders.py"

rule prep_ml_labels:
    input:
        "data/sqlite/activity_data.db"
    output:
        "data/sqlite/activity_data_labelled.db"
    script:
        "scripts/prep_ml_labels.py"

rule train_forest_classifier:
    input:
        "data/sqlite/activity_data_labelled.db"
    output:
        "outputs/forest_classifier_validation.pdf",
        "models/forest_classifier.joblib"
    script:
        "scripts/train_forest_classifier.py"

rule test_forest_classifier:
    input:
        "models/forest_classifier.joblib",
        "data/sqlite/activity_data_labelled.db"
    output:
        "outputs/test_classification_plot.pdf"
    script:
        "scripts/test_forest_classifier.py"

rule train_forest_regressor:
    input:
        "data/sqlite/activity_data_labelled.db"
    output:
        "outputs/forest_regressor_validation.pdf",
        "models/forest_regressor.joblib"
    script:
        "scripts/train_forest_regressor.py"

rule test_forest_regressor:
    input:
        "models/forest_regressor.joblib",
        "data/sqlite/activity_data_labelled.db"
    output:
        "outputs/test_regression_plot.pdf"
    script:
        "scripts/test_forest_regressor.py"

rule plot_regression:
    input:
        "models/linear_model.joblib",
        "data/cleaned/airquality_cleaned.csv"
    output:
        "outputs/reg_plot.pdf"
    script:
        "scripts/plot_linear_model.py"
