#!/usr/bin/env Rscript

# Env ---------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidymodels))
suppressPackageStartupMessages(library(themis))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(Boruta))

# functions ---------------------------------------------------------------

source("workflow/scripts/classification/ml_no_tune.R")

for (i in list.files("workflow/scripts/recipeselectors", full.names = TRUE)) source(i)

# load data ---------------------------------------------------------------

training_data <- read_tsv(snakemake@input[["train"]]) %>% 
  rename(target = all_of(snakemake@params[["target_var"]][["name"]])) %>%
  mutate(target = factor(target)) %>% 
  mutate(target = relevel(target, ref = snakemake@params[["target_var"]][["ref_level"]])) %>%
  rename_with(.fn = function(x) "id", .cols = matches("(^|.*_)id$"))

test_data <- read_tsv(snakemake@input[["test"]]) %>%
  rename(target = all_of(snakemake@params[["target_var"]][["name"]])) %>%
  mutate(target = factor(target)) %>% 
  mutate(target = relevel(target, ref = snakemake@params[["target_var"]][["ref_level"]])) %>%
  rename_with(.fn = function(x) "id", .cols = matches("(^|.*_)id$"))

var_info <- NULL

# set variables -----------------------------------------------------------

ml_alg <- snakemake@wildcards[["ml_alg"]]

model_name <- paste(ml_alg, snakemake@wildcards[["dataset"]], sep = "_")


if (!snakemake@wildcards[["dataset"]] %in% names(snakemake@params[["dataset_proc_steps"]])) {
  steps <- c(snakemake@params[["alg_params"]][["steps"]],
             snakemake@params[["dataset_proc_steps"]][["default"]]) 
} else {
    steps <- c(snakemake@params[["alg_params"]][["steps"]],
               snakemake@params[["dataset_proc_steps"]][[snakemake@wildcards[["dataset"]]]])  
}

if ("n_features" %in% names(snakemake@wildcards)) {
  
  n_features <- as.integer(snakemake@wildcards[["n_features"]])
  steps <- c(steps,
             list("step_select_forests" = list("top_p" = n_features)))
  model_name <- paste(model_name, n_features, sep = "-")
}

print(steps)


# ml ----------------------------------------------------------------------

final_model_res <- 
  run_ml(dat = training_data, 
         target_var = snakemake@params[["target_var"]][["name"]],
         ref_level_target_var = snakemake@params[["target_var"]][["ref_level"]],
         steps = steps,
         algo = snakemake@wildcards[["ml_alg"]],
         test_dat = test_data,
         var_info = var_info)

# get_probs ---------------------------------------------------------------

probs <- final_model_res[["predictions"]]

print("probs gathered")

# get feature importance --------------------------------------------------

if (ml_alg %in% "rf") {
  feature_importances <- extract_fit_parsnip(final_model_res[[1]])$fit$variable.importance %>%
    enframe(name = "feature", value = "feature_importance")
}
if (ml_alg %in% "lasso") {
  feature_importances <- extract_fit_parsnip(final_model_res[[1]]) %>%
    tidy() %>%
    rename(feature = term,
           feature_importance = estimate) 
}
if (ml_alg %in% c("svm", "nnet", "xgb")) {
  feature_importances <- tibble("feature" = NA,
                                "feature_importance" = NA)
}

print("importance gathered")



# get_model_specs ---------------------------------------------------------

if (snakemake@wildcards[["ml_alg"]] %in% c("rf")) {
  model_specs <- 
    extract_fit_parsnip(final_model_res[[1]])$fit[c("num.trees", "num.independent.variables", "mtry", "min.node.size", "prediction.error", "num.samples")] %>%
    bind_cols() %>%
    mutate("preprocessing" = final_model_res[[1]] %>%
             extract_recipe() %>%
             tidy() %>%
             pull(type) %>%
             paste(collapse = ","))
    
} else {
  model_specs <- tibble("num.trees" = NA, 
                        "num.independent.variables" = NA,
                        "mtry" = NA,
                        "min.node.size" = NA, 
                        "prediction.error" = NA,
                        "num.samples" = NA,
                        "reprocessing" = NA)
}

print("specs gathered")



# write -------------------------------------------------------------------

final_model_res[["model"]] %>% 
  write_rds(snakemake@output[["model"]])
probs %>% 
  write_tsv(snakemake@output[["probs"]])
feature_importances %>%
  write_tsv(snakemake@output[["feature_importance"]])
model_specs %>%
  write_tsv(snakemake@output[["model_specs"]])
