#!/usr/bin/env Rscript

# Env ---------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidymodels))
suppressPackageStartupMessages(library(themis))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(Boruta))


if (length(snakemake@log) > 0) {
  sink(snakemake@log[[1]])
}


# functions ---------------------------------------------------------------

source("workflow/scripts/classification/ml_no_tune.R")

for (i in list.files("workflow/scripts/recipeselectors", full.names = TRUE)) source(i)
source("workflow/scripts/classification/filter_prevalence.R")
source("workflow/scripts/classification/filter_lin_comb.R")


# load data ---------------------------------------------------------------

training_data <- read_tsv(snakemake@input[["training"]], col_types = cols()) %>% 
  rename(target = all_of(snakemake@params[["target_var"]][["name"]])) %>%
  mutate(target = factor(target)) %>% 
  mutate(target = relevel(target, ref = snakemake@params[["target_var"]][["ref_level"]])) %>%
  rename_with(.fn = function(x) "id", .cols = matches("(^|.*_)id$"))

if ("var_info" %in% names(snakemake@params)) {
  if (file.exists(snakemake@params[["var_info"]])) {
    var_info <- read_tsv(snakemake@params[["var_info"]])
  } else {
    var_info <- NULL
  }
} else {
  var_info <- NULL
}

if ("add_var" %in% names(snakemake@wildcards)) {
  add_vars_training <- 
    read_tsv(snakemake@input[["add_vars_train"]], col_types = cols())
  
  training_data <-
    training_data %>% 
    rename_with(.fn = function(x) ifelse(x == snakemake@params[["pred"]], "pred", x)) %>% 
    mutate(rep = str_extract(iteration, "^Repeat[:digit:]*"),
           fold = str_extract(iteration, "Fold[:digit:]*$")) %>% 
    select(id, pred, target, rep, fold)
  
  add_vars_assessment <- 
    read_tsv(snakemake@input[["add_vars_assess"]], col_types = cols())
  
  set.seed(8897)
  assessment_data <- 
    read_tsv(snakemake@input[["assessment"]], col_types = cols()) %>% 
    rename_with(.fn = function(x) ifelse(x == snakemake@params[["pred"]], "pred", x)) %>%
    rename_with(.fn = function(x) "id", .cols = matches("(^|.*_)id$")) %>% 
    mutate(rep = str_extract(iteration, "^Repeat[:digit:]*"),
           fold = str_extract(iteration, "Fold[:digit:]*$")) %>% 
    select(id, pred, rep, fold, model) %>% 
    inner_join(add_vars_assessment, by = "id") %>% 
    left_join(training_data %>% select(id, rep, fold) %>% mutate(training_set = TRUE), by = c("id", "rep", "fold")) %>% 
    group_by(id, rep) %>% 
    mutate(training_cat = any(!is.na(training_set))) %>% 
    mutate(assessment_set = sample(5, replace = FALSE)+1) %>% 
    mutate(assessment_cat = (training_set %in% 1 & training_cat) | (assessment_set %in% 2 & !training_cat)) %>% 
    ungroup() %>% 
    filter(assessment_cat) %>% 
    select(-c(training_set, training_cat, assessment_set, assessment_cat))
  
  training_data <-
    training_data %>% 
    inner_join(add_vars_training, by = "id")
  
} else {
  
  assessment_data <- read_tsv(snakemake@input[["assessment"]], col_types = cols()) %>%
    rename(target = all_of(snakemake@params[["target_var"]][["name"]])) %>%
    mutate(target = factor(target)) %>% 
    mutate(target = relevel(target, ref = snakemake@params[["target_var"]][["ref_level"]])) %>%
    rename_with(.fn = function(x) "id", .cols = matches("(^|.*_)id$"))
  
}

# set variables -----------------------------------------------------------



if ("ml_alg" %in% names(snakemake@wildcards)) {
  ml_alg <- snakemake@wildcards[["ml_alg"]]
} else {
  ml_alg <- snakemake@wildcards[["sec_alg"]]
}


model_name <- paste(ml_alg, snakemake@wildcards[["dataset"]], sep = "_")

secondary_model <- "add_var" %in% names(snakemake@wildcards)

sec_no_m_compatible <- TRUE

## Secondary model might be incompatible with non-microbial model
if (secondary_model) {
  if (ncol(add_vars_assessment) <= 2) {
    sec_no_m_compatible <- FALSE
  }
}

if (secondary_model) {
  model_name <- paste(model_name, snakemake@wildcards[["add_var"]], sep = "_")
}

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


if ("fraction" %in% names(snakemake@wildcards)) {
  
  
  fraction <- seq(snakemake@params[["min_sample_frac"]], 
                  1, 
                  length.out = snakemake@params[["sample_fracs"]])[ as.integer(snakemake@wildcards[["fraction"]])]
  
  n_samples <- ceiling(fraction*nrow(training_data))
  
  set.seed(round(fraction*sum(utf8ToInt(snakemake@wildcards[["fraction"]])) + 
                   sum(utf8ToInt(ml_alg)) +
                   sum(utf8ToInt(snakemake@wildcards[["rep"]]))))
  
  training_data <- 
    training_data %>% 
    group_by(target) %>% 
    slice_sample(prop = fraction) %>% 
    ungroup()
  
  model_name <- paste(model_name, snakemake@wildcards[["fraction"]], snakemake@wildcards[["rep"]], sep = "-")
}


# bootstrap ---------------------------------------------------------------

set.seed(8989)
if (snakemake@params[["resampling"]] == "bootstrap") {
  resampled_data <- bootstraps(data = training_data, 
                               times = snakemake@params[["bootstraps"]], 
                               strata = target)
}


# v fold cv ---------------------------------------------------------------

set.seed(8989)
if (snakemake@params[["resampling"]] == "vfold_cv") {
  resampled_data <- vfold_cv(data = training_data, 
                             v = 5, 
                             repeats = 20) 
}



if (snakemake@params[["resampling"]] == "predefined") {
  resampled_data <-
    training_data %>% 
    nest_by(rep, fold, .key = "splits")
  
} else {
  resampled_data <- resampled_data %>% 
    mutate(id = select(., starts_with("id")) %>% apply(1, function(x) paste0(x, collapse = "_")))  
}


# ml ----------------------------------------------------------------------

mods_probs <- lapply(seq(length(resampled_data$splits)), function(i) {
  
  if ("data" %in% names(resampled_data) | "data" %in% names(resampled_data$splits[[i]])) {
    ## cv or bootstrap
    tmp_train <- 
      resampled_data$splits[[i]][["data"]] %>% 
      slice(resampled_data$splits[[i]][["in_id"]])
    tmp_test <- 
      resampled_data$splits[[i]][["data"]] %>% 
      slice(-resampled_data$splits[[i]][["in_id"]])
  } else {
    ## predefined splits
    tmp_train <- resampled_data %>% 
      filter(rep %in% resampled_data$rep[i],
             fold != resampled_data$fold[i]) %>% 
      unnest(splits) %>% 
      ungroup() %>% 
      select(-c(rep, fold))
    tmp_test <- resampled_data$splits[[i]]
  }
  
  tmp_mod_res <- run_ml(dat = tmp_train, 
         ref_level_target_var = snakemake@params[["target_var"]][["ref_level"]],
         iter = as.integer(i), 
         threads = 1,
         steps = steps,
         algo = ml_alg,
         env = NULL,
         test_dat = tmp_test,
         var_info = var_info)
  
  if (secondary_model & sec_no_m_compatible) {
    tmp_mod_res_no_m <-
      run_ml(dat = tmp_train %>% select(-pred), 
             ref_level_target_var = snakemake@params[["target_var"]][["ref_level"]],
             iter = as.integer(i), 
             threads = 1,
             steps = steps,
             algo = ml_alg,
             env = NULL,
             test_dat = tmp_test %>% select(-pred),
             var_info = var_info)
    names(tmp_mod_res_no_m) <- c("no_m_model", "no_m_predictions")
    c(tmp_mod_res, tmp_mod_res_no_m)
  } else {
    tmp_mod_res
  }
})

print("models run")

# get_probs ---------------------------------------------------------------

if (secondary_model) {
  resampled_data <-
    resampled_data %>% 
    mutate(id = paste(rep, fold, sep = "_"))
}

probs <- lapply(seq(length(mods_probs)), function(i) {
  tmp <- mods_probs[[i]][[2]] %>% 
    mutate(iteration = resampled_data$id[i])
  if (secondary_model & sec_no_m_compatible) {
    tmp %>% 
      mutate(no_m = FALSE) %>% 
      bind_rows(mods_probs[[i]][[4]] %>% 
                  mutate(iteration = resampled_data$id[i],
                         no_m = TRUE))
  } else {
    tmp
  }
}) %>% bind_rows() %>% 
  mutate(model = model_name)

print("probs gathered")

# get_assessment_probs ----------------------------------------------------

ass_probs <- lapply(seq(length(mods_probs)), function(i) {
  tmp_n <- names(mods_probs[[i]])
  tmp_n <- tmp_n[ grepl("model", tmp_n)]
  
  lapply(tmp_n, function(mod) {
    
    x <- mods_probs[[i]][[mod]]
    # x <- ml_res[[i]][[2]]  
    
    mod_fit <- extract_fit_parsnip(x)
    
    if (secondary_model) {
      ass <- 
        assessment_data %>% 
        filter(rep %in% resampled_data$rep[i],
               fold %in% resampled_data$fold[i])
    } else {
      ass <- 
        assessment_data
    }
    
    ass <- 
      ass %>% 
      select(which(names(.) %in% extract_recipe(x)$var_info$variable))
    
    train_samples <- extract_recipe(x) %>% 
      bake(new_data = ass) %>% 
      select(id) %>% 
      mutate(training_set = !is.na(id)) %>% 
      mutate(id = ass$id)
    
    if (ml_alg %in% c("xgb", "lasso")) {
      mod_fit %>% 
        predict(new_data = extract_recipe(x) %>% 
                  bake(new_data = ass) %>% select(-any_of(c("id", "target"))), type = "prob") %>% 
        bind_cols(train_samples) %>% 
        mutate(iteration = resampled_data$id[i]) %>% 
        mutate(no_m = mod == "no_m_model")
      
    } else {
      extract_recipe(x) %>% 
        bake(new_data = ass) %>% 
        predict(mod_fit, new_data = ., type = "prob") %>% 
        bind_cols(train_samples) %>% 
        mutate(iteration = resampled_data$id[i]) %>% 
        mutate(no_m = mod == "no_m_model")
    }
    
  }) %>% 
    bind_rows()
  
}) %>% bind_rows() %>% 
  mutate(model = model_name)

print("assessment assessed")

# get feature importance --------------------------------------------------

feature_importances <- lapply(seq(length(mods_probs)), function(i) {
  
  tmp_n <- names(mods_probs[[i]])
  tmp_n <- tmp_n[ grepl("model", tmp_n)]
  
  lapply(tmp_n, function(mod) {
    if (ml_alg %in% "rf") {
      x <- extract_fit_parsnip(mods_probs[[i]][[mod]])$fit$variable.importance %>%
        enframe(name = "feature", value = "feature_importance") %>%
        mutate(iteration = resampled_data$id[i]) %>% 
        mutate(no_m = mod == "no_m_model")
    }
    if (ml_alg %in% c("lasso", "log_reg")) {
      x <- extract_fit_parsnip(mods_probs[[i]][[mod]]) %>%
        tidy() %>%
        rename(feature = term,
               feature_importance = estimate) %>%
        mutate(iteration = resampled_data$id[i]) %>% 
        mutate(no_m = mod == "no_m_model")
    }
    if (ml_alg %in% c("svm", "nnet", "xgb")) {
      x <- tibble("feature" = NA,
                  "feature_importance" = NA,
                  "iteration" = NA)
    }
    x
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows() %>%
  mutate(model = model_name)



# get_model_specs ---------------------------------------------------------

if (ml_alg %in% c("rf")) {
  model_specs <- lapply(seq(length(mods_probs)), function(i) {
    if (ml_alg %in% "rf") {
      x <- extract_fit_parsnip(mods_probs[[i]][[1]])$fit[c("num.trees", "num.independent.variables", "mtry", "min.node.size", "prediction.error", "num.samples")] %>%
        bind_cols() %>%
        mutate("preprocessing" = mods_probs[[i]][[1]] %>%
                 extract_recipe() %>%
                 tidy() %>%
                 pull(type) %>%
                 paste(collapse = ",")) %>%
        mutate(iteration = resampled_data$id[i])
    }
    x
  }) %>% bind_rows() %>%
    mutate(model = model_name)
} else {
  model_specs <- tibble("num.trees" = NA, 
                        "num.independent.variables" = NA,
                        "mtry" = NA,
                        "min.node.size" = NA, 
                        "prediction.error" = NA,
                        "num.samples" = NA,
                        "reprocessing" = NA,
                        "iteration" = NA,
                        "model" = model_name)
}




# write -------------------------------------------------------------------

probs %>% 
  write_tsv(snakemake@output[["probs"]])
print("probs written")
if (model_name %in% "lasso_mags-lididem-fit-2-18") {
  ass_probs %>% 
    slice(0) %>% 
    write_tsv(snakemake@output[["assessment_probs"]])
} else {
  ass_probs %>% 
    write_tsv(snakemake@output[["assessment_probs"]])
}

print("assessment written")
feature_importances %>%
  write_tsv(snakemake@output[["feature_importance"]])
print("importance written")
model_specs %>%
  write_tsv(snakemake@output[["model_specs"]])
print("specs written")

if (length(snakemake@log) > 0) {
  sink()
}
