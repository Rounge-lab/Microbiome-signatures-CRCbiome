

env_vars <- ls()

# prepare data ------------------------------------------------------------

## Lifestyle, diet, demography

lididem <- 
  meta_dat %>% 
  left_join(sample_data %>% select(sample_id, deltaker_id), by = "deltaker_id") %>% 
  select(sample_id, everything(), -deltaker_id) %>% 
  select(sample_id, all_of(variables %>% filter(lididem_crc) %>% pull(var_id))) %>% 
  filter(!is.na(Energi_kcal))

# write function ----------------------------------------------------------

write_training_test <- function(dataset, 
                                dataset_path, 
                                subset_var = "target", 
                                subset_value = c("negative", "positive"), 
                                remove_vars = c("final_result", "train_test", "sel_var", "deltaker_id"),
                                var_info = NULL,
                                create_dummies = FALSE,
                                target_var = "target") {
  dpath <- paste("data/ml-datasets/", dataset_path, "/", sep = "")
  if (!dir.exists(dpath)) dir.create(dpath, recursive = TRUE)
  
  if (create_dummies) {
    dataset <- 
      dataset %>% 
      pivot_longer(cols = where(is.factor) 
                   & !any_of("sample_id", 
                             target_var, 
                             subset_var, 
                             remove_vars, 
                             train_test),
                   names_to = "dummy_names",
                   values_to = "dummy_levels") %>% 
      mutate(dummy_value = 1) %>% 
      pivot_wider(names_from = c(dummy_names, dummy_levels),
                  values_from = dummy_value,
                  values_fill = 0)
      
  }
  
  dataset %>% 
    mutate(new_var = .data[[subset_var]]) %>% 
    filter(new_var %in% subset_value) %>% 
    select(-new_var) %>% 
    filter(train_test %in% "train") %>% 
    select(-any_of(remove_vars)) %>% 
    write_tsv(paste(dpath, "train.tsv", sep = ""))
  
  dataset %>% 
    filter(train_test %in% "train") %>% 
    select(-any_of(remove_vars)) %>% 
    write_tsv(paste(dpath, "assessment_set.tsv", sep = ""))
  
  dataset %>% 
    filter(train_test %in% "test") %>% 
    select(-any_of(remove_vars)) %>% 
    write_tsv(paste(dpath, "test.tsv", sep = ""))
  
  if (!is.null(var_info)) {
    dataset %>% 
      select(-any_of(remove_vars)) %>% 
      names() %>% 
      enframe(value = "var_id") %>% 
      select(-name) %>% 
      left_join(var_info) %>% 
      mutate(var_cat = case_when(is.na(var_cat) ~ "other",
                                 TRUE ~ var_cat)) %>% 
      write_tsv(paste0(dpath, "var_info.tsv"))
  }
  
}


process_dataset <- function(dataset,
                            dataset_name,
                            s_dat = sample_data, 
                            scr_dat = screening_data, 
                            t_t = train_test) {
  dataset %>% 
    left_join(s_dat %>% select(sample_id, deltaker_id), by = "sample_id") %>% 
    left_join(scr_dat %>% select(deltaker_id, detect_worthy_lesions), by = "deltaker_id") %>% 
    left_join(t_t %>% select(deltaker_id, train_test), by = "deltaker_id") %>% 
    select(-deltaker_id) %>% 
    rename(target = detect_worthy_lesions) %>%
    write_training_test(dataset_path = dataset_name)
}

# write datasets ----------------------------------------------------------

## mags
MAGs_abundance %>% 
  process_dataset("mags")

## mphl - (mp4 sgb)
mphlan_abundance %>%
  process_dataset("mphl")

## hum-kegg
ko_hum %>% 
  process_dataset("hum-kegg")
  
## lididem
lididem %>% 
  process_dataset("lididem")
  
## mphl-hum-kegg
mphlan_abundance  %>% 
  left_join(ko_hum) %>% 
  process_dataset("mphl-hum-kegg")

## mags-lididem
MAGs_abundance %>% 
  inner_join(lididem, by = "sample_id") %>% 
  process_dataset("mags-lididem")

## mphl-lididem
mphlan_abundance  %>% 
  inner_join(lididem, by = "sample_id") %>% 
  process_dataset("mphl-lididem")

## hum-kegg-lididem
ko_hum %>% 
  inner_join(lididem, by = "sample_id") %>% 
  process_dataset("hum-kegg-lididem")

## mphl-hum-kegg-lididem
mphlan_abundance  %>% 
  left_join(ko_hum) %>% 
  inner_join(lididem, by = "sample_id") %>% 
  process_dataset("mphl-hum-kegg-lididem")

## Possible to include FIT in model, but this has been omitted here.


# stratifications ---------------------------------------------------------

outcome_splits <- list("neg-v-naa-aa-serr-crc" = str_subset(levels(screening_data$final_result), "^(1|5b|5c|5d|6)."),
                       "neg-nnl-v-aa-serr-crc" = str_subset(levels(screening_data$final_result), "^(1|3|5c|5d|6)."),
                       "neg-v-aa-serr-crc" = str_subset(levels(screening_data$final_result), "^(1|5c|5d|6)."),
                       "neg-v-aa-crc" = str_subset(levels(screening_data$final_result), "^(1|5c|6)."),
                       "neg-v-crc" = str_subset(levels(screening_data$final_result), "^(1|6)."),
                       "neg-v-aa" = str_subset(levels(screening_data$final_result), "^(1|5c)."),
                       "neg-v-serr" = str_subset(levels(screening_data$final_result), "^(1|5d)."))


for(x in c("mphl", "mags", "hum-kegg")) {
  if (x %in% "mphl") {
    tmp <- mphlan_abundance
  }
  if (x %in% "mags") {
    tmp <- MAGs_abundance
  }
  if (x %in% "hum-kegg") {
    tmp <- ko_hum
  }
  
  tmp <- tmp %>% 
    left_join(sample_data %>% select(sample_id, deltaker_id), by = "sample_id") %>% 
    left_join(screening_data %>% select(deltaker_id, detect_worthy_lesions), by = "deltaker_id") %>% 
    left_join(train_test %>% select(deltaker_id, train_test), by = "deltaker_id") %>% 
    rename(target = detect_worthy_lesions)
  
  ## Split by outcome
  for(i in seq(length(outcome_splits))) {
    tmp %>%
      left_join(screening_data %>% select(deltaker_id, final_result), by = "deltaker_id") %>%
      select(-deltaker_id) %>% 
      write_training_test(dataset_path = paste(x, names(outcome_splits)[i], sep = "-"), subset_var = "final_result", subset_value = outcome_splits[[i]])
  }
  ## Split by FIT
  tmp %>%
    left_join(sample_data %>% mutate(sel_var = (FIT_value/5) >= 30) %>% select(sample_id, sel_var), by = "sample_id") %>%
    write_training_test(dataset_path = paste(x, "FIT-high", sep = "-"), subset_var = "sel_var", subset_value = TRUE)
  tmp %>%
    left_join(sample_data %>% mutate(sel_var = (FIT_value/5) >= 30) %>% select(sample_id, sel_var), by = "sample_id") %>%
    write_training_test(dataset_path = paste(x, "FIT-low", sep = "-"), subset_var = "sel_var", subset_value = FALSE)

  ## Split by WCRF
  tmp %>%
    left_join(screening_data %>% select(deltaker_id, id), by = "deltaker_id") %>%
    left_join(meta_dat %>% mutate(sel_var = wcrf_index_main >= 3.5) %>% select(deltaker_id, sel_var), by = "deltaker_id") %>%
    select(-c(deltaker_id)) %>%
    write_training_test(dataset_path = paste(x, "WCRF-high", sep = "-"), subset_var = "sel_var", subset_value = TRUE)
  tmp %>%
    left_join(screening_data %>% select(deltaker_id, id), by = "deltaker_id") %>%
    left_join(meta_dat %>% mutate(sel_var = wcrf_index_main >= 3.5) %>% select(deltaker_id, sel_var), by = "deltaker_id") %>%
    select(-c(deltaker_id)) %>%
    write_training_test(dataset_path = paste(x, "WCRF-low", sep = "-"), subset_var = "sel_var", subset_value = FALSE)

  ## Split by sex
  tmp %>%
    left_join(screening_data %>% mutate(sel_var = kjonn %in% "Male") %>% select(deltaker_id, sel_var), by = "deltaker_id") %>%
    select(-deltaker_id) %>% 
    write_training_test(dataset_path = paste(x, "men", sep = "-"), subset_var = "sel_var", subset_value = TRUE)
  tmp %>%
    left_join(screening_data %>% mutate(sel_var = kjonn %in% "Male") %>% select(deltaker_id, sel_var), by = "deltaker_id") %>%
    select(-deltaker_id) %>% 
    write_training_test(dataset_path = paste(x, "women", sep = "-"), subset_var = "sel_var", subset_value = FALSE)

  ## Split by age
  tmp %>%
    left_join(screening_data %>% mutate(sel_var = age_invitation >= 67) %>% select(deltaker_id, sel_var), by = "deltaker_id") %>%
    select(-deltaker_id) %>% 
    write_training_test(dataset_path = paste(x, "old", sep = "-"), subset_var = "sel_var", subset_value = TRUE)
  tmp %>%
    left_join(screening_data %>% mutate(sel_var = age_invitation >= 67) %>% select(deltaker_id, sel_var), by = "deltaker_id") %>%
    select(-deltaker_id) %>% 
    write_training_test(dataset_path = paste(x, "young", sep = "-"), subset_var = "sel_var", subset_value = FALSE)
  
  ## Split by localization
  tmp %>% 
    left_join(screening_data %>% mutate(sel_var = (distal_acn %in% 1 & proximal_acn %in% 0) | ((distal_acn + proximal_acn) %in% 0)) %>% select(deltaker_id, sel_var), by = "deltaker_id") %>% 
    select(-deltaker_id) %>% 
    write_training_test(dataset_path = paste(x, "distal", sep = "-"), subset_var = "sel_var", subset_value = TRUE)
  tmp %>% 
    left_join(screening_data %>% mutate(sel_var = (distal_acn %in% 0 & proximal_acn %in% 1) | ((distal_acn + proximal_acn) %in% 0)) %>% select(deltaker_id, sel_var), by = "deltaker_id") %>% 
    select(-deltaker_id) %>% 
    write_training_test(dataset_path = paste(x, "proximal", sep = "-"), subset_var = "sel_var", subset_value = TRUE)
  
  ## Split by Center
  tmp %>%
    left_join(screening_data %>% mutate(sel_var = senter %in% "Bærum") %>% select(deltaker_id, sel_var), by = "deltaker_id") %>%
    select(-deltaker_id) %>% 
    write_training_test(dataset_path = paste(x, "Bærum", sep = "-"), subset_var = "sel_var", subset_value = TRUE)
  tmp %>%
    left_join(screening_data %>% mutate(sel_var = senter %in% "Bærum") %>% select(deltaker_id, sel_var), by = "deltaker_id") %>%
    select(-deltaker_id) %>% 
    write_training_test(dataset_path = paste(x, "Moss", sep = "-"), subset_var = "sel_var", subset_value = FALSE)
}


# set up secondary logreg datasets ----------------------------------------

meta_dat_tmp <- 
  meta_dat %>% 
  left_join(sample_data %>% select(sample_id, deltaker_id), by = "deltaker_id") %>% 
  left_join(train_test %>% select(deltaker_id, train_test), by = "deltaker_id") %>% 
  mutate(restr_train_set = final_result %in% c("1. Negative", 
                                               "3. Non neoplastic findings", 
                                               "5c. Advanced adenoma", 
                                               "6. Cancer") & 
           train_test %in% "train",
         FIT_cat = cut(fobt_verdi, breaks = c(15, 20, 25, 35, 75, Inf), 
                       labels = paste0("FIT_", c("15_20", "20_25", "25_35", "35_75", "75_")))) %>% 
  select(sample_id, everything())


fit_dataset <- 
  meta_dat_tmp %>% 
  mutate(logFIT = log(fobt_verdi)) %>% 
  select(sample_id, logFIT, restr_train_set, train_test) %>% 
  rename(id = sample_id)

fit_demo_dataset <- 
  meta_dat_tmp %>% 
  mutate(logFIT = log(fobt_verdi)) %>% 
  select(sample_id, logFIT, kjonn, age_invitation, restr_train_set, train_test) %>% 
  rename(id = sample_id) %>% 
  rename_with(.fn = function(x) variables$var_name[ match(x, variables$var_id)], .cols = which(names(.) %in% variables$var_id))

fit_demo_wcrf_dataset <- 
  meta_dat_tmp %>% 
  mutate(logFIT = log(fobt_verdi)) %>% 
  select(sample_id, logFIT, kjonn, age_invitation, wcrf_index_main, restr_train_set, train_test) %>% 
  rename(id = sample_id) %>% 
  filter(!is.na(wcrf_index_main)) %>% 
  rename_with(.fn = function(x) variables$var_name[ match(x, variables$var_id)], .cols = which(names(.) %in% variables$var_id))

fit_lididem_dataset <-
  meta_dat_tmp %>% 
  mutate(logFIT = log(fobt_verdi)) %>% 
  select(sample_id, logFIT, Utdanning, variables$var_id[variables$lididem_crc], 
         restr_train_set, train_test) %>% 
  rename(id = sample_id) %>% 
  filter(!is.na(Energi_kcal), !is.na(Smoking), !is.na(BMI)) %>% 
  rename_with(.fn = function(x) variables$var_name[ match(x, variables$var_id)], .cols = which(names(.) %in% variables$var_id))

catfit_dataset <-
  meta_dat_tmp %>% 
  select(sample_id, FIT_cat, restr_train_set, train_test) %>% 
  rename(id = sample_id) 

catfit_demo_dataset <- 
  meta_dat_tmp %>% 
  select(sample_id, FIT_cat, kjonn, age_invitation, restr_train_set, train_test) %>% 
  rename(id = sample_id) %>% 
  rename_with(.fn = function(x) variables$var_name[ match(x, variables$var_id)], .cols = which(names(.) %in% variables$var_id))

catfit_demo_wcrf_dataset <- 
  meta_dat_tmp %>% 
  select(sample_id, FIT_cat, kjonn, age_invitation, wcrf_index_main, restr_train_set, train_test) %>% 
  rename(id = sample_id) %>% 
  filter(!is.na(wcrf_index_main)) %>% 
  rename_with(.fn = function(x) variables$var_name[ match(x, variables$var_id)], .cols = which(names(.) %in% variables$var_id))

catfit_lididem_dataset <-
  meta_dat_tmp %>% 
  select(sample_id, FIT_cat, Utdanning, variables$var_id[variables$lididem_crc], 
         restr_train_set, train_test) %>% 
  rename(id = sample_id) %>% 
  filter(!is.na(Energi_kcal), !is.na(Smoking), !is.na(BMI)) %>% 
  rename_with(.fn = function(x) variables$var_name[ match(x, variables$var_id)], .cols = which(names(.) %in% variables$var_id))

restr_PPI <-
  meta_dat_tmp %>% 
  select(sample_id, PPI_antacids_reg_quest_comb, age_invitation, kjonn, FIT_cat,
         restr_train_set, train_test) %>% 
  rename(id = sample_id) %>% 
  rename_with(.fn = function(x) variables$var_name[ match(x, variables$var_id)] %>% 
                str_replace_all(" ", "_") %>% 
                str_replace_all("-", "_"), 
              .cols = which(names(.) %in% variables$var_id))

full_PPI <-
  meta_dat_tmp %>% 
  select(sample_id, PPI_antacids_reg_quest_comb, age_invitation, kjonn, FIT_cat,
         BMI, Sivilstatus_cat2, Alko, symptoms_GI_related_cat, Tarmkreft_Familie,
         restr_train_set, train_test) %>% 
  filter(!is.na(BMI)) %>% 
  rename(id = sample_id) %>% 
  rename_with(.fn = function(x) variables$var_name[ match(x, variables$var_id)] %>% 
                str_replace_all(" ", "_") %>% 
                str_replace_all("-", "_"), 
              .cols = which(names(.) %in% variables$var_id))
  

# write function ----------------------------------------------------------

write_training_test <- function(dataset, 
                                dataset_path, 
                                subset_var = "target", 
                                subset_value = c("negative", "positive"), 
                                remove_vars = c("final_result", "train_test")) {
  dpath <- paste("data/ml-datasets/", dataset_path, "/", sep = "")
  if (!dir.exists(dpath)) dir.create(dpath, recursive = TRUE)
  
  dataset %>% 
    mutate(new_var = .data[[subset_var]]) %>% 
    filter(new_var %in% subset_value) %>% 
    select(-new_var) %>% 
    filter(train_test %in% "train") %>% 
    select(-any_of(remove_vars)) %>% 
    write_tsv(paste(dpath, "train.tsv", sep = ""))
  
  dataset %>% 
    filter(train_test %in% "train") %>% 
    select(-any_of(remove_vars)) %>% 
    write_tsv(paste(dpath, "assessment_set.tsv", sep = ""))
  
  dataset %>% 
    filter(train_test %in% "test") %>% 
    select(-any_of(remove_vars)) %>% 
    write_tsv(paste(dpath, "test.tsv", sep = ""))
  
}

# write datasets ----------------------------------------------------------

## 1. fit_dataset
fit_dataset %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/fit", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))

## 2. fit_demo_dataset
fit_demo_dataset %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/fit_demo", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))

## 3. fit_demo_wcrf_dataset
fit_demo_wcrf_dataset %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/fit_demo_wcrf", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))

## 4. fit_lididem_dataset
fit_lididem_dataset %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/fit_lididem", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))

## 5. catfit_dataset
catfit_dataset %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/catfit", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))

## 6. catfit_demo_dataset
catfit_demo_dataset %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/catfit_demo", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))

## 7. catfit_demo_wcrf_dataset
catfit_demo_wcrf_dataset %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/catfit_demo_wcrf", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))

## 8. catfit_lididem_dataset
catfit_lididem_dataset %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/catfit_lididem", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))

## 9. restr_PPI
restr_PPI %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/restr_PPI", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))

## 10. full_PPI
full_PPI %>% 
  write_training_test(dataset_path = "secondary_lr/norm_nn_v_aa_crc/full_PPI", 
                      subset_var = "restr_train_set", 
                      subset_value = TRUE, 
                      remove_vars = c("restr_train_set", "train_test"))


rm(list = ls()[ !ls() %in% env_vars])
