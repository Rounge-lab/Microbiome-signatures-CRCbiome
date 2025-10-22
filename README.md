# Microbiome-signatures-CRCbiome
This repository contains code for analysis of microbial profiles used to generate results in [Birkeland et al.][medrxiv]

The presented results are generated using R code, which is found in `scripts/analyses/`. Synthetic datasets with individual-level data and files with summary data (such as statistical test outputs) are available in `data/`. The script `scripts/analyses.R` contains the necessary code to recreate figures and tables that are not otherwise available. To run this code, you need to have R and the packages that are listed at the top of the script installed. 

Some calls to perform analyses is commented out, but scripts are available. This includes code to perform preprocessing of individual-level data and to conduct statistical analyses. While the code to conduct statistical testing is functional, running it will overwrite provided the results with associations from synthetic data. 

Code for conducting classification is provided as a snakemake workflow, organized in the snakefile `workflow/classification.smk`. To run this, it is necessary to prep datasets using the `prep_ml_datasets.R` script, which is found in `scripts/analyses/`. Running the snakemake pipeline will carry out the analyses defined in `data/models.tsv`, requires installation of snakemake, with package handling using conda (the software requirements for running the R-encoded scripts are defined in `workflow/envs/tidymodels.yaml`).

The code has been run using Linux CentOS on an HPC system, with sufficient resources to carry out the analyses in a reasonable amount of time. 


[medrxiv]: https://www.medrxiv.org/content/10.1101/2025.10.06.25336873v1