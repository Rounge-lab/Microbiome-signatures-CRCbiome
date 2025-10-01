

# env ---------------------------------------------------------------------

setwd("")

library(tidyverse)
library(pROC)
library(broom)
library(nnet)
library(vegan)
library(ape)
library(gt)
library(gtsummary)
library(ggpubr)
library(cowplot)
library(ggstats)
library(paletteer)


# set variables -----------------------------------------------------------

use_synthetic_data <- TRUE

# functions ---------------------------------------------------------------

source("scripts/analyses/utils.R")

# data preprocessing ------------------------------------------------------

# source("scripts/analyses/preprocess.R")
# preprocess_input_data()

# load data ---------------------------------------------------------------

load_data(synthetic_data = use_synthetic_data)
# source("scripts/analyses/create_synthetic_datasets.R")
# source("scripts/analyses/prep_ml_datasets.R")

# analyses ----------------------------------------------------------------

## Evaluate metaphlan vs MAGs
# source("scripts/analyses/MAG_mphl_comparison.R")

## Conduct alpha and beta diversity testing
# source("scripts/analyses/alpha_beta_host_vars_associations.R")

## Classification DMM groups
# source("scripts/DMM_generation.R")

## Evaluate association btw DMM groups and outcome
# source("scripts/analyses/dmm_outcome_models.R")

## Evaluate ML models
# source("scripts/evaluate_models.R")

## Run host characteristics differential abundance tests adjusted
# source("scripts/analyses/maaslin_meta_vars.R")

## Conduct differential abundance analyses
# source("scripts/analyses/differential_abundance_MaAsLin.R")
# differential_abundance_main_outcome(models = c("main", "strat", "pairwise",
#                                                "low_prev","unadjusted","loc", 
#                                                "colibactin_strat"))

## Evaluate e coli/pks association with outcome
# source("scripts/analyses/ecoli_pks_outcome_eval.R")


# summarize results -------------------------------------------------------

## Summarize participant characteristics
source("scripts/analyses/demo_lifestyle_summary.R")
generate_demo_lifestyle_tables()

## Summarize associations
source("scripts/analyses/plots.R")
plot_microbiome_overview()
plot_loc_associations()
plot_main_da_loc_assoc()
plot_microbial_outcome_associations()
plot_pairwise_beta_div_final_res()
plot_pairwise_diff_abund()
summarize_dmm()
summarize_diff_abund()
plot_strat_da()
plot_ecoli_lit_rep_da()
algorithm_and_dataset_comparison()
summarize_ml_res()

