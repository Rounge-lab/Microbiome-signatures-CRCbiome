
## Functions for model evaluations

evaluate_ml_by_outcome_group <- function(dataset, 
                                         pred_vars = c("primary_prob", "secondary_prob"),
                                         outcome_var = "target",
                                         measures = c("auroc", "sens_spec"),
                                         thresholds = NULL,
                                         strat_outcome = FALSE,
                                         fraction_def_var = NULL) {
  
  
  if (strat_outcome) strats <- c("all", "neg-v-naa", "neg-v-as", "neg-v-aa", "neg-v-crc") else strats <- "all"
  res_tmp <-
    lapply(strats, function(res_subset) {
      if (res_subset %in% "all") tmp <- dataset 
      if (res_subset %in% "neg-v-naa") tmp <- dataset %>% filter(!grepl("^5[cd]", final_result), !grepl("^6", final_result))
      if (res_subset %in% "neg-v-as") tmp <- dataset %>% filter(!grepl("^5[bc]", final_result), !grepl("^6", final_result))
      if (res_subset %in% "neg-v-aa") tmp <- dataset %>% filter(!grepl("^5[bd]", final_result), !grepl("^6", final_result))
      if (res_subset %in% "neg-v-crc") tmp <- dataset %>% filter(!grepl("^5", final_result))
      
      tmp %>% 
        mutate(selection = res_subset) %>% 
        get_stat(pred_vars = pred_vars,
                 outcome_var = outcome_var,
                 thresholds = thresholds,
                 measure = measures,
                 fraction_def_var = fraction_def_var)
      
        
    }) %>% 
    bind_rows()
  
  return(res_tmp)
}

define_threshold <- function(dataset,
                             pred_vars,
                             def_sens_var_id = "final_result",
                             def_sens_var_lvl = "6. Cancer",
                             def_sens_val = 0.9) {
  
  dataset %>% 
    rename(sel_var = all_of(def_sens_var_id)) %>% 
    pivot_longer(all_of(pred_vars), names_to = "predictor", values_to = "prediction") %>%
    group_by(rep, fold, predictor) %>% 
    filter(str_detect(sel_var, def_sens_var_lvl)) %>% 
    summarize(threshold = quantile(prediction, 1-def_sens_val), 
              n_in_group = n(),
              .groups = "drop")
}

get_stat <- function(dataset, 
                     pred_vars = c("primary_prob", "secondary_prob", "fobt_verdi"),
                     outcome_var = "target",
                     measure = c("auroc", "sens_spec", "fractions"),
                     thresholds = NULL,
                     fraction_def_var = NULL) {
  
  if ("sens_spec" %in% measure & is.null(thresholds)) {
    stop("Need to specify thresholds when evaluating sensitivity and specificity")
  }
  if ("n_by_outcome" %in% measure & (is.null(fraction_def_var) | is.null(thresholds))) {
    stop("Need to specify thresholds and which variable to enumerate
         number of positive/negative predictions when requesting 'n_by_outcome'
         as an outcome measure.")
  }
  if ("n_by_outcome" %in% measure) {
    if (!fraction_def_var %in% names(dataset)) {
      stop("Variable to enumerate by needs to be available.")
    }
  }
  
  if (!is.null(fraction_def_var)) {
    if (outcome_var == fraction_def_var) {
      dataset <-
        dataset %>% 
        bind_cols(dataset %>% 
                    select(fraction_def_var = all_of(fraction_def_var)))
      fraction_def_var <- "fraction_def_var"
    }  
  }
  
  dataset <- 
    dataset %>% 
    rename(target = all_of(outcome_var)) %>% 
    pivot_longer(all_of(pred_vars), names_to = "predictor", values_to = "prediction")
  
  lapply(measure, function(measure_i) {
    if (measure_i == "auroc") {
      dataset %>% 
        group_by(predictor, rep, fold) %>% 
        summarize(auroc = ifelse(all(c("negative", "positive") %in% target) &
                                   !all(is.na(prediction)), 
                                 roc(target~prediction, levels = c("negative", "positive"), direction = "<")$auc[[1]],
                                 NA),
                  .groups = "drop")  
    } else if (measure_i == "sens_spec") {
      dataset %>% 
        left_join(thresholds, by = join_by(predictor, rep, fold)) %>% 
        mutate(prediction = prediction > threshold) %>% 
        group_by(predictor, rep, fold) %>% 
        summarize(sensitivity = sum(prediction & target %in% "positive")/sum(target %in% "positive"),
                  specificity = sum(!prediction & target %in% "negative")/sum(target %in% "negative"),
                  .groups = "drop")  
    } else if (measure_i == "n_by_outcome") {
      dataset %>% 
        left_join(thresholds, by = join_by(predictor, rep, fold)) %>% 
        rename(enum_var = all_of(fraction_def_var)) %>% 
        mutate(prediction = prediction > threshold) %>% 
        group_by(predictor, rep, fold) %>% 
        nest() %>% 
        mutate(enum = map(.x = data, .f = function(dat = .x) {
          dat %>% 
            count(prediction, enum_var)
        })) %>% 
        select(-data)
    }
  }) %>% 
    reduce(left_join, by = c("predictor", "rep", "fold")) %>% 
    left_join(dataset %>% 
                ungroup() %>% 
                select(predictor, rep, fold, any_of("selection")) %>% 
                distinct(),
              by = join_by(predictor, rep, fold))
}

get_stat_by_strat <- function(dataset, 
                              ..., strat_var) {
  dataset %>% 
    rename(strat_var = all_of(strat_var)) %>% 
    nest_by(strat_var) %>% 
    mutate(preds = list(get_stat(data, ...))) %>% 
    select(-data) %>% 
    ungroup() %>% 
    unnest(preds) %>% 
    rename(strat_val = strat_var) %>% 
    mutate(strat_var = strat_var)
}

evaluate_models <- function(metadata, dataset_name, primary_algorithm, 
                            secondary_algorithm, secondary_add_var, 
                            included_in_training = "all",
                            primary_eval = FALSE, strat_outcome = FALSE,
                            outcome_cat = "detect_worthy_lesions",
                            meta_for_threshold = NULL,
                            measures = c("auroc", "sens_spec"),
                            fraction_def_var = NULL,
                            def_threshold_by_var = "final_result",
                            def_threshold_by = "6. Cancer",
                            def_by_threshold = 0.9) {
  
  tmp_secondary_probs <-
    read_tsv(paste0("data/models/", primary_algorithm, "_", dataset_name, 
                    "/training_sec/norm_nn_v_aa_crc/", secondary_add_var, 
                    "/", secondary_algorithm, "_probs.tsv"), col_types = cols())
  
  tmp_secondary_assessment_probs <-
    read_tsv(paste0("data/models/", primary_algorithm, "_", dataset_name, 
                    "/training_sec/norm_nn_v_aa_crc/", secondary_add_var, 
                    "/", secondary_algorithm, "_assessment_probs.tsv"), col_types = cols()) %>% 
    select(-training_set) %>% 
    left_join(tmp_secondary_probs %>% 
                select(id, iteration) %>% 
                distinct() %>% 
                mutate(secondary_training_set = TRUE), 
              by = join_by(id, iteration)) %>% 
    mutate(secondary_training_set = ifelse(is.na(secondary_training_set), FALSE, TRUE)) %>% 
    mutate(rep = str_extract(iteration, "Repeat[:digit:]*"),
           fold = str_extract(iteration, "Fold[:digit:]*")) %>% 
    rename(secondary_prob = .pred_positive) %>% 
    rename(secondary_model = model) %>% 
    mutate(primary_model = paste0(primary_algorithm, "_", dataset_name)) %>% 
    select(-.pred_negative) %>% 
    (function(x) {
      if ("no_m" %in% names(x)) {
        x %>% 
          mutate(no_m = ifelse(no_m, "only_add_var", "secondary_prob")) %>% 
          pivot_wider(names_from = no_m, values_from = secondary_prob)
      } else {
        x
      }
    })
  
  if (primary_eval) {
    
    ## Primary model evaluation
    
    tmp_primary_probs <-
      read_tsv(paste0("data/models/", primary_algorithm, "_", dataset_name, 
                      "/training/probs.tsv"), col_types = cols()) %>% 
      rename(primary_model = model) %>% 
      select(id, iteration, primary_model) %>% 
      mutate(primary_training_set = TRUE)
    
    tmp_primary_assessment_probs <-
      read_tsv(paste0("data/models/", primary_algorithm, "_", dataset_name, 
                      "/training/assessment_probs.tsv"), col_types = cols()) %>% 
      rename(primary_model = model) %>% 
      inner_join(tmp_secondary_assessment_probs %>% 
                   select(id, iteration, primary_model) %>% 
                   distinct(), 
                 by = join_by(id, iteration, primary_model)) %>% 
      select(-training_set) %>% 
      left_join(tmp_primary_probs, 
                by = join_by(id, iteration, primary_model)) %>% 
      mutate(primary_training_set = ifelse(is.na(primary_training_set), FALSE, TRUE)) %>% 
      mutate(rep = str_extract(iteration, "Repeat[:digit:]*"),
             fold = str_extract(iteration, "Fold[:digit:]*")) %>% 
      rename(primary_prob = .pred_positive)
    
    prim_train_set <- included_in_training
    
    if (!(prim_train_set == "external" & 
          all(tmp_primary_assessment_probs$primary_training_set))) {
      
      tmp_res <- 
        tmp_primary_assessment_probs %>% 
        mutate(train_set_filter = case_when(prim_train_set == "all" ~ TRUE,
                                            prim_train_set == "internal" ~ primary_training_set,
                                            prim_train_set == "external" ~ !primary_training_set)) %>% 
        filter(train_set_filter) 
      
      if (is.null(meta_for_threshold)) {
        tmp_meta <- metadata
      } else {
        tmp_meta <- meta_for_threshold
      }
      
      threshold_df <-
        tmp_res %>% 
        inner_join(tmp_meta, by = "id") %>% 
        define_threshold(pred_vars = c("primary_prob"),
                         def_sens_var_id = def_threshold_by_var,
                         def_sens_var_lvl = def_threshold_by,
                         def_sens_val = def_by_threshold)
      
      tmp_res <- 
        tmp_res %>% 
        inner_join(metadata, by = "id")
      
        
      tmp_res %>% 
        reframe(evaluate_ml_by_outcome_group(dataset = ., 
                                             pred_vars = c("primary_prob"), 
                                             outcome_var = outcome_cat, 
                                             thresholds = threshold_df, 
                                             strat_outcome = strat_outcome,
                                             measures = measures,
                                             fraction_def_var = fraction_def_var)) %>% 
        mutate(included_in_training = included_in_training) %>% 
        left_join(tmp_res %>% select(primary_model, rep, fold) %>% distinct(), by = join_by(rep, fold))
      
    }
  } else {
    
    ## Secondary model evaluation
    sec_train_set <- included_in_training
    if (!(sec_train_set == "external" & 
          all(tmp_secondary_assessment_probs$secondary_training_set))) {
      
      tmp_res <- 
        tmp_secondary_assessment_probs %>% 
        mutate(train_set_filter = case_when(sec_train_set == "all" ~ TRUE,
                                            sec_train_set == "internal" ~ secondary_training_set,
                                            sec_train_set == "external" ~ !secondary_training_set)) %>% 
        filter(train_set_filter) 
      
      if (is.null(meta_for_threshold)) {
        tmp_meta <- metadata
      } else {
        tmp_meta <- meta_for_threshold
      }
      
      p_vars <- c("secondary_prob", "FIT_value", "only_add_var")
      p_vars <- p_vars[ p_vars %in% c(names(tmp_res), names(metadata), names(tmp_meta))]
      
      threshold_df <-
        tmp_res %>% 
        inner_join(tmp_meta, by = "id") %>% 
        define_threshold(pred_vars = p_vars,
                         def_sens_var_id = def_threshold_by_var,
                         def_sens_var_lvl = def_threshold_by,
                         def_sens_val = def_by_threshold)
      
      tmp_res <-
        tmp_res %>% 
        inner_join(metadata, by = "id")
      
      tmp_res %>% 
        reframe(evaluate_ml_by_outcome_group(dataset = ., 
                                             pred_vars = p_vars, 
                                             outcome_var = outcome_cat, 
                                             thresholds = threshold_df, 
                                             strat_outcome = strat_outcome,
                                             measures = measures,
                                             fraction_def_var = fraction_def_var)) %>% 
        mutate(included_in_training = included_in_training) %>% 
        left_join(tmp_res %>% select(primary_model, secondary_model, rep, fold) %>% distinct(), by = join_by(rep, fold))
    }
  }
}

combine_prim_sec <- function(pri_res, sec_res) {
  pri_res %>% 
    mutate(secondary_model = sec_res$secondary_model[1]) %>% 
    select(-included_in_training) %>% 
    bind_rows(sec_res %>% 
                select(-included_in_training) %>% 
                mutate(predictor = case_when(FALSE ~ "only_add_var",
                                             TRUE ~ predictor)))
}

get_test_set_evaluation <- function(output = "AUROC") {
  
  final_models <- list("ml" = c("rf"),
                       "dataset" = c("mphl", "mags"),
                       "add_var" = c("fit", "fit_demo_wcrf"))
  
  hold_out_locs <-
    sapply(final_models$ml, function(ml) {
      sapply(final_models$dataset, function(ds) {
        sapply(final_models$add_var, function(av) {
          paste0("data/models/final_model/", ml, "_", 
                 ds, "/secondary_lr/norm_nn_v_aa_crc/", av,"_probs.tsv")
        })
      })
    }) %>% 
    as.vector()
  
  test_set_classification <- 
    hold_out_locs %>% 
    lapply(function(hold_out_loc) {
      read_tsv(hold_out_loc, col_types = cols())
    }) %>% 
    bind_rows()
  
  # final model -------------------------------------------------------------
  if (output == "roc") {
    
    d <-
      test_set_classification %>% 
      rename(sample_id = id) %>% 
      left_join(sample_data %>% 
                  left_join(screening_data, by = "deltaker_id") %>% 
                  select(sample_id, detect_worthy_lesions, FIT_value) %>% 
                  rename(target = detect_worthy_lesions), 
                by = "sample_id") %>% 
      pivot_longer(c(primary_pred, secondary_pred, FIT_value), names_to = "predictor") %>% 
      filter(!(predictor %in% "primary_pred" & no_m),
             !(str_detect(predictor, "FIT") & no_m)) %>% 
      group_by(predictor, primary_model, add_var, no_m, secondary_alg) %>%
      nest() %>% 
      mutate(roc_obj = map(.x = data, .f = function(x = .x) {
        x %>% 
          roc_plot_obj(response = "target", predictor = "value")
      })) %>% 
      select(-data) %>% 
      unnest(roc_obj)
  }
  
  if (output == "FIT_sens_spec") {
    fit_cat <-
      sample_data %>% 
      left_join(screening_data %>% select(deltaker_id, detect_worthy_lesions), by = "deltaker_id") %>% 
      select(sample_id, detect_worthy_lesions, FIT_value) %>% 
      rename(target = detect_worthy_lesions) %>% 
      mutate(FIT_over_20 = (FIT_value/5)>20,
             FIT_over_25 = (FIT_value/5)>25,
             FIT_over_35 = (FIT_value/5)>35,
             FIT_over_75 = (FIT_value/5)>75) %>% 
      mutate(across(starts_with("FIT"), .fns = function(x) x*1)) %>% 
      select(-FIT_value)
    
    ## Sensitivities and specificities associated with particular FIT cutoffs  
    d <- 
      test_set_classification %>% 
      filter(!no_m) %>% 
      rename(sample_id = id) %>% 
      select(sample_id, primary_model, add_var, secondary_alg) %>% 
      left_join(fit_cat, by = "sample_id") %>% 
      pivot_longer(starts_with("FIT"), names_to = "predictor") %>% 
      group_by(predictor, primary_model, add_var, secondary_alg) %>%
      nest() %>% 
      mutate(roc_obj = map(.x = data, .f = function(x = .x) {
        x %>% 
          roc_plot_obj(response = "target", predictor = "value")
      })) %>% 
      select(-data) %>% 
      unnest(roc_obj) %>% 
      filter(sensitivities != 1,
             specificities != 1,
             !is.na(thresholds))
  }
  
  if (output == "AUROC") {
    
    d <-
      test_set_classification %>% 
      rename(sample_id = id) %>% 
      left_join(sample_data %>% 
                  left_join(screening_data, by = "deltaker_id") %>% 
                  select(sample_id, detect_worthy_lesions, FIT_value) %>% 
                  rename(target = detect_worthy_lesions), 
                by = "sample_id") %>%
      pivot_longer(c(primary_pred, secondary_pred, FIT_value), names_to = "predictor") %>% 
      filter(!(predictor %in% "primary_pred" & no_m),
             !(str_detect(predictor, "FIT") & no_m)) %>%
      group_by(predictor, primary_model, add_var, no_m, secondary_alg) %>%
      nest() %>% 
      mutate(roc_obj = map(.x = data, .f = function(x = .x) {
        tmp <- 
          roc(x$target~x$value, levels = c("negative", "positive"), direction = "<") %>% 
          ci.auc(method = "bootstrap") %>% 
          as.vector()
        tibble("auroc" = tmp[2],
               "ci_low" = tmp[1],
               "ci_high" = tmp[3])
      })) %>% 
      select(-data) %>% 
      unnest(roc_obj)
  }
  
  if (output == "sens_spec") {
    
    test_set_probs <-
      test_set_classification %>% 
      rename(sample_id = id) %>% 
      left_join(sample_data %>% 
                  left_join(screening_data, by = "deltaker_id") %>% 
                  select(sample_id, detect_worthy_lesions, FIT_value) %>% 
                  rename(target = detect_worthy_lesions), 
                by = "sample_id") %>%
      pivot_longer(c(primary_pred, secondary_pred, FIT_value), names_to = "predictor") %>% 
      filter(!(predictor %in% "primary_pred" & no_m),
             !(str_detect(predictor, "FIT") & no_m))
    
    define_threshold_by <- "all_75"
    # define_threshold_by <- "CRC_90"
    
    training_thresholds <- 
      lapply(c("mphl", "mags"), function(prim_dataset) {
        lapply(c("fit", "fit_demo_wcrf"), function(add_var) {
            read_tsv(paste0("data/models/rf_", 
                            prim_dataset, 
                            "/training_sec/norm_nn_v_aa_crc/",
                            add_var, "/log_reg_assessment_probs.tsv"), 
                     col_types = cols()) %>%
            select(-c(.pred_negative, training_set)) %>%
            (function(x) {
              if ("no_m" %in% names(x)) {
                x %>% 
                  mutate(no_m = ifelse(no_m, "only_add_var", "secondary_pred")) %>% 
                  pivot_wider(names_from = no_m, values_from = .pred_positive)
              } else {
                x %>% 
                  rename(secondary_pred = .pred_positive)
              }
            }) %>% 
            left_join(read_tsv(paste0("data/models/rf_", 
                                      prim_dataset, 
                                      "/training/assessment_probs.tsv"), 
                               col_types = cols()) %>% 
                        select(primary_pred = .pred_positive, id, iteration, primary_model = model),
                      by = join_by(id, iteration)) %>% 
            left_join(sample_data %>% 
                        select(id = sample_id, FIT_value) %>% 
                        mutate(FIT_value = FIT_value/5), 
                      by = join_by(id)) %>% 
            pivot_longer(any_of(c("primary_pred", "secondary_pred", "FIT_value", "only_add_var")),
                         names_to = "predictor", values_to = "prediction") %>% 
            (function(x) {
              if (define_threshold_by == "CRC_90") {
                x %>% 
                  inner_join(sample_data %>% 
                               select(id = sample_id, deltaker_id) %>% 
                               inner_join(screening_data %>% 
                                            select(deltaker_id, final_result) %>% 
                                            filter(str_detect(final_result, "6. Cancer")), 
                                          by = join_by(deltaker_id)) %>% 
                               select(-deltaker_id),
                             by = join_by(id)) %>% 
                  group_by(predictor, model, primary_model, iteration) %>% 
                  summarize(threshold = quantile(prediction, 0.1), .groups = "drop") %>% 
                  group_by(predictor, model, primary_model) %>% 
                  summarize(mean_threshold = mean(threshold), .groups = "drop")
              } else if (define_threshold_by == "all_75") {
                x %>% 
                  group_by(predictor, model, primary_model, iteration) %>% 
                  summarize(threshold = quantile(prediction, 0.25), .groups = "drop") %>% 
                  group_by(predictor, model, primary_model) %>% 
                  summarize(mean_threshold = mean(threshold), .groups = "drop")
              } 
            })
            
        }) %>% 
          bind_rows()
      }) %>% 
      bind_rows() %>% 
      mutate(add_var = str_remove(model, "log_reg_(mphl|mags)_"))
    
    d <-
      test_set_probs %>%
      mutate(value = case_when(predictor %in% "FIT_value" ~ value/5,
                               TRUE ~ value)) %>% 
      mutate(predictor = ifelse(no_m, "only_add_var", predictor)) %>% 
      select(-no_m) %>% 
      group_by(primary_model, add_var, secondary_alg, predictor) %>% 
      nest() %>% 
      left_join(training_thresholds %>% select(-model)) %>% 
      mutate(ci_sens_spec = map2(.x = data, .y = mean_threshold, .f = function(preds = .x, thr = .y) {
        roc_obj <- pROC::roc(preds$target ~ preds$value, levels = c("negative", "positive"), direction = "<")
        tmp <- pROC::ci.thresholds(roc_obj, thresholds = thr)
        
        tmp$specificity %>% 
          as.data.frame() %>% 
          pivot_longer(everything(), 
                       names_to = "CI", 
                       values_to = "specificity") %>% 
          left_join(tmp$sensitivity %>% 
                      as.data.frame() %>% 
                      pivot_longer(everything(), 
                                   names_to = "CI", 
                                   values_to = "sensitivity"), by = join_by(CI)) %>% 
          mutate(CI = case_when(CI == "2.5%" ~ "CI_low",
                                CI == "97.5%" ~ "CI_high",
                                TRUE ~ "bootstrap_median")) %>% 
          pivot_longer(c(sensitivity, specificity), names_to = "measure", values_to = "value") %>% 
          pivot_wider(names_from = CI, values_from = value)
      })) %>% 
      select(-data) %>% 
      unnest(ci_sens_spec)
    
  }
  
  return(d)
}

