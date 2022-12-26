This repository contains datasets, model code and notebooks used for all experiments in the [PROSTATA: Protein Stability Assessment using
Transformers](https://www.biorxiv.org/content/10.1101/2022.12.25.521875v1) paper.

Folders:

- DATA contains the datasets used for PROSTATA training and testing in the format used by the datasets authors. Also the dataset introduced in this article is available [here](https://github.com/mitiau/PROSTATA/blob/main/DATA/dataset_our_w_clusters_v2.5.pkl) .
- DATASETS contains the same datasets converted to a format used for model training.
- PDB contains the PDB files downloaded during conversion.

Notebooks:

- 00.generate_datasets.ipynb - Process the DATA. folder and generate the DATASETS folder
- 01.test_models_by_folds.ipynb - Test each individual model in the ensemble using5-fold cross validation
- 02.test_models_on_other_datasets_ensemble.ipynb - test the PROSTATA ensemble on various combinations of train and test datasets. Output is logged to the PROSTATA_experiments_pearson.log
- 03.train_final_ensemble.ipynb - train the ensemble on all data for the online tool.
- PROSTATA_tool.ipynb - Colab notebook for PROSTATA. Predict DDG Values for single mutation on user sequence.

Other files:

- PROSTATA_experiments_pearson.log - Logs of experiments run by 02.test_models_on_other_datasets_ensemble.ipynb notebook.
- LICENSE - Apache License 2.0
- Readme.MD - This file
