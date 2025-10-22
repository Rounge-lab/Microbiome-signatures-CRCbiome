
# env ---------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# load data ---------------------------------------------------------------

# source("workflow/scripts/analyses/load_data.R")
source("scripts/analyses/load_synthetic_data.R")

## Secondary model specs
fi_training <- read_tsv(paste0("data/models/",
                               snakemake@wildcards[["ml_alg"]], 
                               "_", 
                               snakemake@wildcards[["dataset"]], 
                               "/training_sec/",
                               snakemake@wildcards[["preselection"]], 
                               "/",
                               snakemake@wildcards[["add_var"]], 
                               "/log_reg_feature_importance.tsv"), col_types = cols())

# functions ---------------------------------------------------------------

create_dummy <- function(values) {
  tmp_lvls <- levels(as.factor(values))
  tmp_dum <- lapply(tmp_lvls, function(lv_i) (values %in% lv_i)*1)
  names(tmp_dum) <- tmp_lvls
  tmp_dum %>% bind_cols()
}

predict_testset <- function(alg, dataset, add_var) {
  
  model_name <- paste(alg, dataset, sep = "_")
  
  ## Secondary model addvars
  add_vars <- read_tsv(paste0("data/ml-datasets/secondary_lr/norm_nn_v_aa_crc/", add_var, "/test.tsv"), col_types = cols())
  
  ## Test set predictions
  testset_preds <- read_tsv(paste0("data/models/final_model/", model_name, "/probs.tsv"), col_types = cols())
  

  ## recode add_vars 
  
  tmp_add_vars <-
    add_vars %>% 
    ## Dummy character variables to match "factor_level" pattern
    mutate(across(.cols = where(is.character) & !matches("id"), .fns = create_dummy, .unpack = TRUE)) %>% 
    select(-(where(is.character) & !matches("id")))
  
  ## restrict dataset
  tmp_fi_training <-
    fi_training %>% 
    filter(!is.na(feature_importance))
  
  if (!"no_m" %in% names(tmp_fi_training)) tmp_fi_training$no_m <- FALSE
  
  lapply(unique(tmp_fi_training$no_m), function(no_m_cat) {
    ## define secondary model
    tmp_model_coefs <-
      tmp_fi_training %>% 
      filter(no_m %in% no_m_cat) %>% 
      select(-c(std.error, statistic, p.value)) %>% 
      group_by(feature) %>% 
      ## Use the mean coefficient from the training set models
      summarize(feature_importance = mean(feature_importance), .groups = "drop") %>% 
      pivot_wider(names_from = feature, values_from = feature_importance)# %>%
    
    intercept <- tmp_model_coefs %>% pull(`(Intercept)`)
    
    tmp_model_coefs <- tmp_model_coefs %>% select(-`(Intercept)`)
    
    tmp_rename_vars = list("for modeling" = c(`pred` = ".pred_positive"),
                           "for output" = c(`primary_pred` = "pred"))
    
    ## Predict
    testset_preds %>% 
      rename(any_of(tmp_rename_vars[["for modeling"]])) %>% 
      inner_join(tmp_add_vars, by = "id") %>% 
      select(id, tmp_model_coefs %>% names()) %>% 
      mutate(new_lincomb = (as.matrix(.[,-1]) %*% as.vector(t(tmp_model_coefs)) + intercept) %>% as.data.frame() %>%  pull(1),
             secondary_pred = 1/(1+exp(-new_lincomb))) %>% 
      select(-new_lincomb) %>% 
      left_join(testset_preds %>% select(id, target), by = "id") %>% 
      select(everything(), secondary_pred) %>% 
      rename(any_of(tmp_rename_vars[["for output"]])) %>% 
      mutate(primary_model = model_name,
             add_var = add_var,
             no_m = no_m_cat,
             secondary_alg = "log_reg")
  }) %>% 
    bind_rows()
}

# predict -----------------------------------------------------------------

secondary_preds <-
  predict_testset(alg = snakemake@wildcards[["ml_alg"]],
                  dataset = snakemake@wildcards[["dataset"]],
                  add_var = snakemake@wildcards[["add_var"]]) %>%
  select(id, primary_pred, secondary_pred, primary_model, add_var, no_m, secondary_alg)

secondary_preds %>% 
  write_tsv(snakemake@output[["secondary_predictions"]])
