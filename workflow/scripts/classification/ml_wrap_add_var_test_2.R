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

training_data <- 
  read_tsv(snakemake@input[["probs"]]) %>%
  rename(target = all_of(snakemake@params[["target"]])) %>%
  rename(pred = all_of(snakemake@params[["pred"]])) %>%
  mutate(target = factor(target)) %>%
  mutate(target = relevel(target, ref = snakemake@params[["ref_level"]])) %>%
  rename_with(.fn = function(x) "id", .cols = matches("(^|.*_)id$")) %>%
  select(id, pred, target, iteration, model) %>% 
  mutate(rep = str_extract(iteration, "Repeat[:digit:]*"),
         fold = str_extract(iteration, "Fold[:digit:]*"))


add_vars_training <- 
  read_tsv(snakemake@input[["add_vars_train"]])


# set variables -----------------------------------------------------------

steps <- snakemake@params[["steps"]]
if (any(names(steps) %in% "step_normalize")) {
  steps <- steps[-which(names(steps) %in% "step_normalize")]  
}

# ml ----------------------------------------------------------------------

cl <- makeCluster(snakemake@threads)

clusterExport(cl = cl, varlist = c("snakemake", "steps", "run_ml"))

add_var_models <-
  training_data %>% 
  inner_join(add_vars_training) %>% 
  group_by(model) %>% 
  group_split() %>% 
  parLapply(cl = cl, X = ., function(x) {
    
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(tidymodels))
    suppressPackageStartupMessages(library(themis))
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(Boruta))
    
    
    tmp_train <- 
      x %>% 
      select(-c(iteration, model, rep, fold))
    
    tmp_test <- 
      x %>% 
      select(-c(iteration, model, rep, fold))
    
    tmp <- run_ml(dat = tmp_train, 
                  target_var = snakemake@params[["target"]],
                  ref_level_target_var = snakemake@params[["ref_level"]],
                  steps = steps,
                  algo = "log_reg",
                  test_dat = tmp_test)
    
    list("model" = tmp[["model"]],
         "probs" = tmp[["predictions"]],
         "model_name" = x$model[1])
  })

stopCluster(cl)

print("models run")


# load assessment data ----------------------------------------------------

assessment_data <- 
  read_tsv(snakemake@input[["assessment_probs"]]) %>%
  rename(pred = all_of(snakemake@params[["pred"]])) %>%
  rename_with(.fn = function(x) "id", .cols = matches("(^|.*_)id$")) %>% 
  select(id, pred, iteration, model)

add_vars_assessment <- 
  read_tsv(snakemake@input[["add_vars_assess"]])


# get_robs ---------------------------------------------------------------

probs <- lapply(seq(length(add_var_models)), function(i) {
  
  tmp_mod <- add_var_models[[i]][["model_name"]]
  tmp_iters <- unique(training_data$iteration[ training_data$model %in% tmp_mod])
  ## Make prediction per iteration (even if the model is independent of iteration)
  ## This is so we retain the iteration information
  lapply(seq(length(tmp_iters)), function(ii) {
    
    x <- add_var_models[[i]][[1]]
    
    mod_fit <- extract_fit_parsnip(x)
    
    dat <- training_data %>% 
      filter(model %in% add_var_models[[i]][["model_name"]],
             iteration %in% tmp_iters[ii]) %>% 
      inner_join(add_vars_assessment, by = "id") %>% 
      select(which(names(.) %in% extract_recipe(x)$var_info$variable))
    
    if (nrow(dat) > 0) {
      extract_recipe(x) %>%
        bake(new_data = dat) %>%
        predict(mod_fit, new_data = ., type = "prob") %>%
        bind_cols(dat %>% select(id)) %>% 
        mutate(iteration = tmp_iters[ii],
               model = add_var_models[[i]][["model_name"]])
    }
    
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows() %>% 
  left_join(training_data %>% rename(original_prob = pred) %>%  select(model, iteration, id, original_prob))


print("probs gathered")

# get_assessment_probs ----------------------------------------------------


## Make resampling of new samples consistent with resampling scheme
reps <- length(unique(training_data$rep))
folds <- length(unique(training_data$fold))

cl <- makeCluster(snakemake@threads)
clusterExport(cl = cl, varlist = c("reps", "folds"))

## Create one fold per sample, "folds" number of folds, repeated "reps" times
assessment_resampling <- 
  assessment_data %>%
  select(id, model) %>% 
  unique() %>% 
  left_join(training_data %>% select(id, model) %>% unique() %>% mutate(dat = "training")) %>% 
  filter(is.na(dat)) %>% 
  group_by(model) %>% 
  group_split() %>% 
  parLapply(cl = cl, X = ., function(x) {
    
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(tidymodels))
    
    tmp_strat <- 
      x %>% 
      vfold_cv(v = folds, repeats = reps)
    
    if (reps > 1) {
      lapply(tmp_strat$splits, function(xx) {
        xx$data %>% 
          select(id) %>% 
          slice(-xx$in_id) %>% 
          mutate(rep = xx$id$id,
                 fold = xx$id$id2,
                 model = x$model[1],
                 iteration = paste(rep, fold, sep = "_"))
      }) %>% 
        bind_rows()  
    } else {
      lapply(tmp_strat$splits, function(xx) {
        xx$data %>% 
          select(id) %>% 
          slice(-xx$in_id) %>% 
          mutate(rep = NA,
                 fold = xx$id$id,
                 model = x$model[1],
                 iteration = fold)
      }) %>% 
        bind_rows()  
    }
  }) %>% 
  bind_rows()

stopCluster(cl)

if (nrow(assessment_resampling) > 0) {
  
  ass_probs <- lapply(seq(length(add_var_models)), function(i) {
    tmp_mod <- add_var_models[[i]][["model_name"]]
    tmp_iters <- unique(training_data$iteration[ training_data$model %in% tmp_mod])
    ## Make prediction per iteration (even if the model is independent of iteration)
    ## This is so we retain the iteration information
    lapply(seq(length(tmp_iters)), function(ii) {
    # lapply(seq(length(add_var_models[[i]])), function(ii) {
      x <- add_var_models[[i]][[1]]
      
      mod_fit <- extract_fit_parsnip(x)
      
      ass <- 
        assessment_data %>% 
        filter(model %in% add_var_models[[i]][["model_name"]]) %>% 
        inner_join(assessment_resampling %>% filter(model %in% add_var_models[[i]][["model_name"]]), by = c("id", "iteration", "model"))  %>% 
        filter(iteration %in% tmp_iters[ii]) %>% 
        inner_join(add_vars_assessment, by = "id") %>% 
        select(which(names(.) %in% extract_recipe(x)$var_info$variable))
      
      if (nrow(ass) > 0) {
        extract_recipe(x) %>%
          bake(new_data = ass) %>%
          predict(mod_fit, new_data = ., type = "prob") %>%
          bind_cols(ass %>% select(id)) %>% 
          mutate(iteration = tmp_iters[ii],
                 model = add_var_models[[i]][["model_name"]])
      }
      
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows() %>% 
    left_join(assessment_data %>% rename(original_prob = pred) %>%  select(model, iteration, id, original_prob))
} else {
  ass_probs <- 
    tibble(".pred_negative" = character(),
           ".pred_positive" = character(),
           "id" = character(),
           "iteration" = character(),
           "model" = character(),
           "original_prob" = character())
}

print("assessment assessed")

# # get feature importance --------------------------------------------------

feature_importances <- 
  lapply(seq(length(add_var_models)), function(i) {
    extract_fit_parsnip(add_var_models[[i]][[1]]) %>%
      tidy() %>%
      rename(feature = term,
             feature_importance = estimate) %>%
      mutate(model = add_var_models[[i]][["model_name"]])
  }) %>% 
  bind_rows()


# write -------------------------------------------------------------------

probs %>% 
  write_tsv(snakemake@output[["probs"]])
ass_probs %>%
  write_tsv(snakemake@output[["assessment_probs"]])
feature_importances %>%
  write_tsv(snakemake@output[["fi"]])
