

run_ml <- function(dat,
                   target_var = "target",
                   ref_level_target_var = "negative",
                   training_prop = 0.8,
                   iter = 1,
                   threads = 1,
                   steps,
                   importance_measure = "permutation",
                   env = NULL,
                   test_dat = NULL,
                   algo,
                   var_info = NULL,
                   verbose = FALSE) {
  
  dat <- dat %>% 
    rename(target = all_of(target_var)) %>% 
    mutate(target = factor(target) %>% relevel(ref = ref_level_target_var)) %>% 
    rename_with(.fn = function(x) "id", .cols = matches("(^|.*_)id$"))
  
  set.seed(iter)
  if (is.null(test_dat)) {
    tr_te_set <-  dat %>% 
      initial_split(prop = training_prop, strata = target) 
  } else {
    test_dat <- test_dat %>% 
      rename(target = all_of(target_var)) %>% 
      mutate(target = factor(target))
    
    combined <- bind_rows(dat, test_dat)
    ind <- list(analysis = seq(nrow(dat)), assessment = nrow(dat) + seq(nrow(test_dat)))
    tr_te_set <- make_splits(ind, combined)
  }
  
  dat_train <- tr_te_set %>% 
    training()
  
  dat_test <- tr_te_set %>% 
    testing()
  
  # Set up model ------------------------------------------------------------
  
  set.seed(iter*52)
  ## Define a recipe:
  rec <- dat_train %>% 
    recipe() %>% 
    update_role(everything(), new_role = "predictor") %>% 
    update_role(all_of("target"), new_role = "outcome")
  
  if (any(grepl("(^|.*_)id$", names(dat)))) {
    rec <- rec %>% 
      update_role(matches("(^|.*_)id$"), new_role = "id")
  }
  
  rec <- rec %>% 
    step_novel(all_nominal_predictors())
  
  for (i in seq(length(steps))) {
    
    if(!is.null(var_info)) {
      if(steps[[i]] %in% var_info$var_cat) {
        sel_vars <- var_info %>% 
          filter(var_cat %in% steps[[i]]) %>% 
          pull(var_id)
      } else {
        sel_vars <- dat %>% 
          select(-grep("(^|.*_)id$", names(.)), - target) %>% 
          names() %>% 
          enframe() %>% 
          pull(value)
      }
    } else {
      sel_vars <- dat %>% 
        select(-grep("(^|.*_)id$", names(.)), - target) %>% 
        names() %>% 
        enframe() %>% 
        pull(value)
    }
    
    if (names(steps)[i] == "step_filter_lin_comb") {
      
      rec <- 
        rec %>% 
        step_select_lin_comb(
          all_predictors(),
          outcome = "target",
          n_top_contr = steps[["step_filter_lin_comb"]][["n_top_contr"]],
          selection_criterion = steps[["step_filter_lin_comb"]][["selection_criterion"]],
          adjustment_dataset = steps[["step_filter_lin_comb"]][["adjustment_dataset"]],
          min_lrt_p = steps[["step_filter_lin_comb"]][["min_lrt_p"]])
    }
    
    if (names(steps)[i] == "step_prevalence") {
      
      max_z <- 1-steps[["step_prevalence"]][["min_p"]]
      min_z <- 1-steps[["step_prevalence"]][["max_p"]]
      
      rec <- rec %>%
        step_filter_zero_fraction(all_numeric_predictors(), 
                                  min_zero_fraction = min_z, 
                                  max_zero_fraction = max_z)
      
    }
        
    
    if (names(steps)[i] == "step_select_boruta") {
      ## Boruta for different groups of variables    
      if(!is.null(var_info)) {
        if(all(steps[[i]] %in% var_info$var_cat)) {
          for (ii in seq(length(steps[[i]]))) {
            sel_vars <- var_info %>% 
              filter(var_cat %in% steps[[i]][[ii]]) %>% 
              pull(var_id)
            rec <- rec %>% 
              step_select_boruta(any_of(!!sel_vars), outcome = "target")
          }
        } else {
          rec <- rec %>% 
            step_select_boruta(all_predictors(), outcome = "target")
        }
      } else {
        rec <- rec %>% 
          step_select_boruta(all_predictors(), outcome = "target")
      }
    }
    
    
    if (names(steps)[i] == "step_smote") {
      rec <- rec %>% 
        step_smote(target)
    }
    
    if (names(steps)[i] == "step_adasyn") {
      
      ## over ratio cannot be lower than the ratio between lowest and highest in the training set
      min_over_ratio <- dat_train %>% 
        count(target) %>% 
        summarize(obs_over_ratio = min(n)/max(n)) %>% 
        mutate(min_over_ratio = min(obs_over_ratio*1.1, 1)) %>% 
        pull(min_over_ratio)
      ## Neighbours must be set low enough to allow small sample sets
      max_neighbours <- dat_train %>% 
        count(target) %>% 
        mutate(max_neighbours = n-1) %>% 
        pull(max_neighbours) %>% 
        min()
      sample_neighbours <- min(5,max_neighbours)
      if (steps[[i]][["over_ratio"]] > min_over_ratio) {
        
        
        if (!(any(names(steps) == "step_dummy"))) {
          print("variables dummified to enable adasyn")
          rec <- rec %>% 
            step_dummy(all_nominal(), -all_outcomes(), -has_role("id"))
        } else {
          if (which(names(steps) == "step_dummy") > i ) {
            print("variable dummification moved to enable adasyn")
            rec <- rec %>% 
              step_dummy(all_nominal(), -all_outcomes(), -has_role("id"))
            
            steps <- steps[ names(steps) != "step_dummy"]
          } 
        }
        
        rec <- rec %>% 
          step_adasyn(target, seed = iter,
                      over_ratio = steps[[i]][["over_ratio"]],
                      neighbors = sample_neighbours)
      } else {
        if (verbose) {
          print("over ratio set too low - will not perform adasyn")  
        }
      }
      
    }
    
    if (names(steps)[i] == "step_relevel") {
      rec <- rec %>% 
        step_relevel(target, 
                     ref_level = ref_level_target_var)
    }
    
    if (names(steps)[i] == "step_nzv") {
      rec <- rec %>% 
        step_nzv(all_predictors(), 
                 unique_cut = steps[[i]][["unique_cut"]],
                 freq_cut = steps[[i]][["freq_cut"]])
    }
    
    if (names(steps)[i] == "step_dummy") {
      rec <- rec %>% 
        step_dummy(all_nominal(), -all_outcomes(), -has_role("id"))
    }
    
    if (names(steps)[i] == "step_select_roc") {
      if ("top_p" %in% names(steps[[i]])) {
        rec <- rec %>% 
          step_select_roc(all_numeric_predictors(), 
                          outcome = "target", 
                          top_p = steps[[i]][["top_p"]])
      } else {
        if ("threshold" %in% names(steps[[i]])) {
          rec <- rec %>% 
            step_select_roc(all_numeric_predictors(),
                            outcome = "target",
                            threshold = steps[[i]][["threshold"]])
        }
      }
    }
    
    if (names(steps)[i] == "step_select_forests") {
      if ("top_p" %in% names(steps[[i]])) {
        rec <- rec %>% 
          step_select_forests(all_predictors(), 
                          outcome = "target", 
                          top_p = steps[[i]][["top_p"]])
      } else {
        if ("threshold" %in% names(steps[[i]])) {
          rec <- rec %>% 
            step_select_forests(all_predictors(),
                            outcome = "target",
                            threshold = steps[[i]][["threshold"]])
        }
      }
    }
    
    if (names(steps)[i] == "step_impute_median_mode") {
      rec <- rec %>% 
        step_impute_median(all_numeric_predictors()) %>% 
        step_impute_mode(all_nominal_predictors())
    }
    
    if (names(steps)[i] == "step_naomit") {
      rec <- rec %>% 
        step_naomit(everything(), skip = TRUE)
    }
    
    if (names(steps)[i] == "step_downsample") {
      rec <- rec %>% 
        step_downsample(target)
    }
    
    if (names(steps)[i] == "step_normalize") {
      rec <- rec %>% 
        step_normalize(all_numeric(), 
                       -all_outcomes(), 
                       -has_role("id"))
    }
    
    if (names(steps)[i] == "step_corr") {
      rec <- rec %>% 
        step_corr(all_numeric(),
                  threshold = steps[[i]][["threshold"]],
                  method = steps[[i]][["method"]])
    }
  }
  
  if (algo == "rf") {
    
    ml_mod <- rand_forest(mode = "classification",
                          mtry = round(sqrt(.cols())), 
                          trees = 1000, 
                          min_n = 10) %>% 
      set_engine("ranger",
                 importance = !!importance_measure,
                 probability = TRUE,
                 seed = iter^2)
    
  }
  
  if (algo == "xgb") {
    ml_mod <- boost_tree(mode = "classification",
                         mtry = round(sqrt(.cols())),
                         trees = 1000, 
                         tree_depth = 4, learn_rate = 0.3,
                         loss_reduction = 0) %>% 
      set_engine("xgboost", 
                 importance = !!importance_measure,
                 probability = TRUE,
                 seed = iter^2)
  }
  
  if (algo == "lasso") {
    ml_mod <- logistic_reg(
      mode = "classification",
      penalty = 1e-3,
      mixture = 1
    ) %>% 
      set_engine("glmnet")
  }
  
  if (algo == "svm") {
    ml_mod <- svm_rbf(
      mode = "classification",
      cost = 2,
      rbf_sigma = 1e-3
    ) %>% 
      set_engine("kernlab", )
  }
  
  if (algo == "log_reg") {
    if (length(unique(dat_train$target)) > 2) {
      ml_mod <- multinom_reg(engine = "nnet",
                             mode = "classification")
    } else {
      ml_mod <- logistic_reg(engine = "glm",
                             mode = "classification")  
    }
  }
  
  if (algo == "nnet") {
    ml_mod <- mlp(mode = "classification",
                  activation = "softmax",
                  engine = "nnet",
                  hidden_units = 1,
                  dropout = 0,
                  epochs = 25)
  }
  
  ml_wf <- workflow() %>% 
    add_recipe(rec) %>% 
    add_model(ml_mod)
  
  ml_fit <- ml_wf %>% 
    fit(data = dat_train)
  
  ml_predict <- ml_fit %>% 
    predict(new_data = dat_test, type = "prob") %>% 
    bind_cols(dat_test %>% select(id, target))
  
  return(list("model" = ml_fit, "predictions" = ml_predict))
}
