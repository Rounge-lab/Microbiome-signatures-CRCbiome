


models <- read_tsv("data/models.tsv", col_types = cols())

tmp_meta <- 
  sample_data %>% 
  select(id = sample_id, deltaker_id) %>% 
  left_join(screening_data %>% 
              select(deltaker_id, final_result, detect_worthy_lesions, final_result_cat_neg, distal_acn, proximal_acn, FIT_value = fobt_verdi),
            by = "deltaker_id")

source("workflow/scripts/classification/model_eval_functions.R")

# model comparisons -------------------------------------------------------

alg_comp_assessments <-
  models %>% 
  filter(model_type %in% "base") %>% 
  left_join(expand.grid(dataset = models$dataset, 
                        sec_add_var = "fit", 
                        sec_alg = "log_reg") %>% 
              tibble(), by = "dataset", relationship = "many-to-many") %>% 
  rowwise() %>% 
  mutate(primary_assessment = list(evaluate_models(metadata = tmp_meta, 
                                                     dataset_name = dataset,
                                                     primary_algorithm = alg,
                                                     secondary_algorithm = sec_alg,
                                                     secondary_add_var = sec_add_var,
                                                     primary_eval = TRUE,
                                                     strat_outcome = FALSE))) 

dir.create("data/models/evaluations/alg_comp_models/", showWarnings = FALSE)

alg_comp_assessments %>% 
  select(dataset, alg, primary_assessment) %>%
  unnest(primary_assessment) %>%
  write_tsv("data/models/evaluations/alg_comp_models/primary_assessments.tsv")

# feature strat training --------------------------------------------------

## Run feature strat models only for microbial and FIT value (continuous)
feature_strat_training_res <- 
  models %>% 
  filter(model_type %in% "feature-strat") %>% 
  left_join(expand.grid(dataset = models$dataset, 
                        sec_add_var = "fit", 
                        sec_alg = c("log_reg")) %>% 
              tibble(), by = "dataset") %>% 
  slice_head(n = 2) %>% 
  rowwise() %>% 
  mutate(secondary_assessment = list(evaluate_models(metadata = tmp_meta, 
                                                     dataset_name = dataset,
                                                     primary_algorithm = alg,
                                                     secondary_algorithm = sec_alg,
                                                     secondary_add_var = sec_add_var,
                                                     strat_outcome = FALSE))) %>% 
  mutate(primary_assessment = list(evaluate_models(metadata = tmp_meta, 
                                                   dataset_name = dataset,
                                                   primary_algorithm = alg,
                                                   secondary_algorithm = sec_alg,
                                                   secondary_add_var = sec_add_var,
                                                   primary_eval = TRUE,
                                                   strat_outcome = FALSE))) %>% 
  mutate(primary_internal = list(evaluate_models(metadata = tmp_meta, 
                                                 dataset_name = dataset,
                                                 primary_algorithm = alg,
                                                 secondary_algorithm = sec_alg,
                                                 secondary_add_var = sec_add_var,
                                                 primary_eval = TRUE,
                                                 included_in_training = "internal",
                                                 strat_outcome = FALSE))) %>% 
  mutate(primary_external = list(evaluate_models(metadata = tmp_meta, 
                                                 dataset_name = dataset,
                                                 primary_algorithm = alg,
                                                 secondary_algorithm = sec_alg,
                                                 secondary_add_var = sec_add_var,
                                                 primary_eval = TRUE,
                                                 included_in_training = "external",
                                                 strat_outcome = FALSE))) %>% 
  mutate(assessments = list(combine_prim_sec(pri_res = primary_assessment, sec_res = secondary_assessment))) %>% 
  select(-c(secondary_assessment, primary_assessment)) %>% 
  ungroup()
  
dir.create("data/models/evaluations/feature_strat_models/", showWarnings = FALSE)

feature_strat_training_res %>% 
  select(-starts_with("primary")) %>% 
  unnest(assessments) %>% 
  write_tsv("data/models/evaluations/feature_strat_models/assessments.tsv")

feature_strat_training_res %>% 
  select(-c(assessments, primary_external)) %>% 
  unnest(primary_internal) %>% 
  write_tsv("data/models/evaluations/feature_strat_models/primary_assessments_internal.tsv")

feature_strat_training_res %>% 
  select(-c(assessments, primary_internal)) %>% 
  unnest(primary_external) %>% 
  write_tsv("data/models/evaluations/feature_strat_models/primary_assessments_external.tsv")


# main models -------------------------------------------------------------

main_models_outcome_strat <-
  models %>% 
  filter(dataset %in% c("mags", "mphl"),
         alg %in% "rf") %>% 
  left_join(expand.grid(dataset = .$dataset, 
                        sec_add_var = c("fit", "fit_demo", "fit_demo_wcrf"), 
                        sec_alg = c("log_reg")) %>% 
              tibble(), by = "dataset", relationship = "many-to-many") %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(secondary_assessment = list(evaluate_models(metadata = tmp_meta,
                                                     dataset_name = dataset,
                                                     primary_algorithm = alg,
                                                     secondary_algorithm = sec_alg,
                                                     secondary_add_var = sec_add_var,
                                                     strat_outcome = TRUE,
                                                     measures = c("auroc", "sens_spec"),
                                                     fraction_def_var = "final_result"))) %>%
  mutate(primary_assessment = list(evaluate_models(metadata = tmp_meta,
                                                   dataset_name = dataset,
                                                   primary_algorithm = alg,
                                                   secondary_algorithm = sec_alg,
                                                   secondary_add_var = sec_add_var,
                                                   primary_eval = TRUE,
                                                   strat_outcome = TRUE,
                                                   measures = c("auroc", "sens_spec"),
                                                   fraction_def_var = "final_result"))) %>%
  ## Do not include exclusively proximal in distal evaluation, and exclude non-advanced lesions
  mutate(distal_sec_assessment = list(evaluate_models(metadata = tmp_meta %>%
                                                        filter(!(distal_acn == 0 & proximal_acn == 1),
                                                               !str_detect(final_result, "^5b.")) %>%
                                                        mutate(distal_acn = ifelse(distal_acn == 1, "positive", "negative")),
                                                  dataset_name = dataset,
                                                  primary_algorithm = alg,
                                                  secondary_algorithm = sec_alg,
                                                  secondary_add_var = sec_add_var,
                                                  primary_eval = FALSE,
                                                  strat_outcome = FALSE,
                                                  outcome_cat = "distal_acn"))) %>%
  ## Do not include exclusively distal in proximal evaluation, and exclude non-advanced lesions
  mutate(proximal_sec_assessment = list(evaluate_models(metadata = tmp_meta %>%
                                                          filter(!(distal_acn == 1 & proximal_acn == 0),
                                                                 !str_detect(final_result, "^5b.")) %>%
                                                          mutate(proximal_acn = ifelse(proximal_acn == 1, "positive", "negative")),
                                                  dataset_name = dataset,
                                                  primary_algorithm = alg,
                                                  secondary_algorithm = sec_alg,
                                                  secondary_add_var = sec_add_var,
                                                  primary_eval = FALSE,
                                                  strat_outcome = FALSE,
                                                  outcome_cat = "proximal_acn"))) %>%
  mutate(distal_prim_assessment = list(evaluate_models(metadata = tmp_meta %>%
                                                         filter(!(distal_acn == 0 & proximal_acn == 1),
                                                                !str_detect(final_result, "^5b.")) %>%
                                                         mutate(distal_acn = ifelse(distal_acn == 1, "positive", "negative")),
                                                       dataset_name = dataset,
                                                       primary_algorithm = alg,
                                                       secondary_algorithm = sec_alg,
                                                       secondary_add_var = sec_add_var,
                                                       primary_eval = TRUE,
                                                       strat_outcome = FALSE,
                                                       outcome_cat = "distal_acn"))) %>%
  mutate(proximal_prim_assessment = list(evaluate_models(metadata = tmp_meta %>%
                                                           filter(!(distal_acn == 1 & proximal_acn == 0),
                                                                  !str_detect(final_result, "^5b.")) %>%
                                                           mutate(proximal_acn = ifelse(proximal_acn == 1, "positive", "negative")),
                                                         dataset_name = dataset,
                                                         primary_algorithm = alg,
                                                         secondary_algorithm = sec_alg,
                                                         secondary_add_var = sec_add_var,
                                                         primary_eval = TRUE,
                                                         strat_outcome = FALSE,
                                                         outcome_cat = "proximal_acn"))) %>%
  mutate(assessments = list(combine_prim_sec(pri_res = primary_assessment, sec_res = secondary_assessment))) %>% 
  mutate(distal_assessments = list(combine_prim_sec(pri_res = distal_prim_assessment, sec_res = distal_sec_assessment))) %>%
  mutate(proximal_assessments = list(combine_prim_sec(pri_res = proximal_prim_assessment, sec_res = proximal_sec_assessment))) %>%
  rename(primary_dataset = dataset, primary_alg = alg) %>% 
  ungroup() %>% 
  select(-ends_with("assessment"))
  
dir.create("data/models/evaluations/main_models/", showWarnings = FALSE)  

main_models_outcome_strat %>% 
  select(primary_alg, primary_dataset, sec_add_var, sec_alg, assessments) %>% 
  unnest(assessments) %>% 
  # select(-enum) %>% 
  write_tsv("data/models/evaluations/main_models/assessments.tsv")

main_models_outcome_strat %>% 
  select(primary_alg, primary_dataset, sec_add_var, sec_alg, distal_assessments) %>% 
  unnest(distal_assessments) %>% 
  write_tsv("data/models/evaluations/main_models/distal_assessments.tsv")
main_models_outcome_strat %>% 
  select(primary_alg, primary_dataset, sec_add_var, sec_alg, proximal_assessments) %>% 
  unnest(proximal_assessments) %>% 
  write_tsv("data/models/evaluations/main_models/proximal_assessments.tsv")


# Enumerate classified cases ----------------------------------------------

evaluate_models(metadata = tmp_meta,
                dataset_name = "mphl",
                primary_algorithm = "rf",
                secondary_algorithm = "log_reg",
                secondary_add_var = "fit_demo_wcrf",
                strat_outcome = FALSE,
                measures = c("n_by_outcome"),
                fraction_def_var = "detect_worthy_lesions",
                def_threshold_by_var = "detect_worthy_lesions",
                def_threshold_by = "(negative|positive)",
                def_by_threshold = 0.8) %>% 
  unnest(enum)


main_models_enumerate_by_threshold <-
  models %>% 
  filter(dataset %in% c("mags", "mphl"),
         alg %in% "rf") %>% 
  left_join(expand.grid(dataset = .$dataset, 
                        sec_add_var = c("fit", "fit_demo", "fit_demo_wcrf"), 
                        sec_alg = c("log_reg")) %>% 
              tibble(), by = "dataset", relationship = "many-to-many") %>% 
  distinct() %>% 
  mutate(a = "a") %>% 
  left_join(tibble(threshold = seq(0.5, 0.9, by = 0.05),
                   # def_outcome = "^[:digit:]",
                   def_outcome = "(Negative|>= 3 Non-advanced adenomas|Advanced serrated|Advanced adenoma|Cancer)",
                   fraction_def_var = "final_result_cat_neg",
                   # fraction_def_var = "final_result",
                   a = "a") %>% 
              bind_rows(tibble(threshold = seq(0.5, 0.9, by = 0.05),
                               def_outcome = "(negative|positive)",
                               fraction_def_var = "detect_worthy_lesions",
                               a = "a")), 
            by = "a", relationship = "many-to-many") %>% 
  select(-a) %>% 
  rowwise() %>% 
  mutate(secondary_assessment = list(evaluate_models(metadata = tmp_meta,
                                                     dataset_name = dataset,
                                                     primary_algorithm = alg,
                                                     secondary_algorithm = sec_alg,
                                                     secondary_add_var = sec_add_var,
                                                     strat_outcome = FALSE,
                                                     measures = c("n_by_outcome"),
                                                     fraction_def_var = fraction_def_var,
                                                     def_threshold_by_var = fraction_def_var,
                                                     def_threshold_by = def_outcome,
                                                     def_by_threshold = threshold))) %>%
  mutate(primary_assessment = list(evaluate_models(metadata = tmp_meta,
                                                   dataset_name = dataset,
                                                   primary_algorithm = alg,
                                                   secondary_algorithm = sec_alg,
                                                   secondary_add_var = sec_add_var,
                                                   primary_eval = TRUE,
                                                   strat_outcome = FALSE,
                                                   measures = c("n_by_outcome"),
                                                   fraction_def_var = fraction_def_var,
                                                   def_threshold_by_var = fraction_def_var,
                                                   def_threshold_by = def_outcome,
                                                   def_by_threshold = threshold))) %>%
  mutate(assessments = list(combine_prim_sec(pri_res = primary_assessment, 
                                                     sec_res = secondary_assessment))) %>% 
  rename(primary_dataset = dataset, primary_alg = alg) %>% 
  ungroup() %>% 
  select(-ends_with("assessment"))

main_models_enumerate_by_threshold %>% 
  mutate(def_outcome = "all") %>% 
  select(primary_alg, primary_dataset, sec_add_var, sec_alg, fraction_def_var, def_outcome, threshold, assessments) %>%
  unnest(assessments) %>% 
  unnest(enum) %>% 
  write_tsv("data/models/evaluations/main_models/outcome_enumerations_tmp.tsv")


# stratified predictions --------------------------------------------------

read_tsv("data/lmr/crcbiome_comorbid_prescriptions.tsv", col_types = cols()) %>% 
  filter(within %in% "12 months",
         drug_cat != "antibiotics_pres") %>% 
  select(deltaker_id = id, drug_cat, prescribed) %>% 
  pivot_wider(names_from = drug_cat, values_from = prescribed)
  

strat_dat <-
  sample_data %>% 
  left_join(load_dmm() %>% select(sample_id, dmm = gr), by = "sample_id") %>% 
  left_join(screening_data %>% select(deltaker_id, senter, kjonn, age_cat, id, detect_worthy_lesions, final_result), by = "deltaker_id") %>% 
  left_join(meta_dat %>% select(deltaker_id, Antacids = PPI_antacids_reg_quest_comb, fobt_verdi,
                                antibiotics_use = antibiotics_reg_quest_comb,
                                colo_hemorrhoids, colo_diverticulitis), by = "deltaker_id") %>% 
  left_join(meta_dat_cat %>% select(deltaker_id, wcrf_index_main, Smoking, Utdanning), by = "deltaker_id") %>% 
  left_join(read_tsv("data/lmr/crcbiome_comorbid_prescriptions.tsv", col_types = cols()) %>% 
              filter(within %in% "12 months",
                     drug_cat != "antibiotics_pres") %>% 
              select(deltaker_id = id, drug_cat, prescribed) %>% 
              pivot_wider(names_from = drug_cat, values_from = prescribed), by = join_by(deltaker_id)) %>% 
  mutate(FIT_value_cat = cut(fobt_verdi, breaks = c(0, 25, 40, 80, Inf), labels = c("15-25 µg/g", "25-35 µg/g", "35-75 µg/g", "75+ µg/g"))) %>%
  mutate(across(ends_with("_pres"), function(x) ifelse(x, "Drug dispensed", "No drug dispensed"))) %>% 
  select(id = sample_id, dmm, senter, kjonn, antibiotics_use, Antacids, FIT_value_cat, FIT_value = fobt_verdi, 
         age_cat, Smoking, Utdanning, WCRF = wcrf_index_main,
         colo_hemorrhoids, colo_diverticulitis,
         Diabetes = diab_2_pres, `Cardiovascular disease` = cvd_pres, `Chronic obstructive pulmonary disease` = copd_pres, 
         detect_worthy_lesions, final_result)

strat_pred_results <-
  models %>% 
  filter(dataset %in% c("mags", "mphl"),
         alg %in% "rf") %>% 
  left_join(expand.grid(dataset = .$dataset, 
                        sec_add_var = c("fit", "fit_demo", "fit_demo_wcrf"), 
                        sec_alg = c("log_reg")) %>% 
              tibble(), by = "dataset", relationship = "many-to-many") %>% 
  mutate(a = "a") %>% 
  left_join(strat_dat %>% 
              select(-c(id, detect_worthy_lesions, final_result, FIT_value)) %>% 
              pivot_longer(everything(), names_to = "strat_var", values_to = "strat_val") %>% 
              distinct() %>% 
              filter(!strat_val %in% c("Missing"), !is.na(strat_val)) %>% 
              mutate(a = "a"), by = "a", relationship = "many-to-many") %>% 
  select(-a) %>%
  rowwise() %>% 
  mutate(tmp_strat_dat = map2(.x = strat_var, .y = strat_val, function(var, val) {
    strat_dat %>% 
      rename(strat_var = all_of(var)) %>% 
      filter(strat_var %in% val) %>% 
      inner_join(train_test %>% 
                   filter(train_test %in% "train") %>% 
                   inner_join(sample_data %>% 
                                select(id = sample_id, deltaker_id), 
                              by = "deltaker_id") %>% 
                   select(id), 
                 by = "id")
  })) %>% 
  mutate(secondary_assessment = list(evaluate_models(metadata = tmp_strat_dat,
                                                     dataset_name = dataset,
                                                     primary_algorithm = alg,
                                                     secondary_algorithm = sec_alg,
                                                     secondary_add_var = sec_add_var,
                                                     primary_eval = FALSE,
                                                     strat_outcome = FALSE,
                                                     meta_for_threshold = tmp_meta,
                                                     measures = c("auroc")))) %>% 
  mutate(primary_assessment = list(evaluate_models(metadata = tmp_strat_dat,
                                                   dataset_name = dataset,
                                                   primary_algorithm = alg,
                                                   secondary_algorithm = sec_alg,
                                                   secondary_add_var = sec_add_var,
                                                   primary_eval = TRUE,
                                                   strat_outcome = FALSE,
                                                   meta_for_threshold = tmp_meta,
                                                   measures = c("auroc")))) %>% 
  mutate(assessments = list(combine_prim_sec(pri_res = primary_assessment, sec_res = secondary_assessment))) %>% 
  ungroup() %>% 
  mutate(n_in_group = map_int(.x = tmp_strat_dat, .f = function(x = .x) nrow(x))) %>% 
  select(-c(tmp_strat_dat, ends_with("assessment")))


dir.create("data/models/evaluations/stratified_predictions/", showWarnings = FALSE)

strat_pred_results %>% 
  unnest(assessments) %>% 
  write_tsv("data/models/evaluations/stratified_predictions/assessments.tsv")


# Evaluate final models ---------------------------------------------------

dir.create("data/models/final_model/evaluations", showWarnings = FALSE)

get_test_set_evaluation(output = "AUROC") %>% 
  write_tsv("data/models/final_model/evaluations/AUROC.tsv")
get_test_set_evaluation(output = "roc") %>% 
  write_tsv("data/models/final_model/evaluations/roc.tsv")
get_test_set_evaluation(output = "FIT_sens_spec") %>% 
  write_tsv("data/models/final_model/evaluations/FIT_sens_spec.tsv")
get_test_set_evaluation(output = "sens_spec") %>% 
  write_tsv("data/models/final_model/evaluations/sens_spec.tsv")

  

