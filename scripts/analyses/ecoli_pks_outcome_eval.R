
# ## Diff abund of pks-ecoli
# colibactin_strat_diff_abund <- 
#   read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_colibactin_strat_models.tsv", 
#            col_types = cols()) %>% 
#   mutate(sign_lab = case_when(qval < 0.05 ~ "   *",
#                               TRUE ~ "")) %>% 
#   filter(str_detect(feature, "ecoli"))

colibactin <- read_tsv("data/input_processed/pks_detection.tsv", col_types = cols())

adj_vars <- 
  sample_data %>% 
  select(sample_id, deltaker_id) %>% 
  left_join(meta_dat %>% 
              select(deltaker_id, senter, 
                     age_invitation, kjonn, 
                     wcrf_index_main, Smoking, 
                     colo_hemorrhoids, colo_IBD,
                     antibiotics_reg_quest_comb,
                     Utdanning, detect_worthy_lesions,
                     final_result_cat_neg), 
            by = "deltaker_id") %>% 
  select(-c(deltaker_id))


cb_for_mod <-
  colibactin %>% 
  select(sample_id, ecoli_pos, pks_classification) %>% 
  left_join(adj_vars) %>%
  select(-sample_id) %>% 
  mutate(Smoking = fct_na_level_to_value(Smoking, "Missing")) %>%
  mutate(Utdanning = fct_na_level_to_value(Smoking, "Missing")) %>%
  # count(Smoking, detect_worthy_lesions)
  mutate(across(where(is.character), as.factor)) %>% 
  mutate(pks = paste(ecoli_pos, pks_classification)) 


evaluate_ecoli_assoc <- function(dat, outcome, pks_strat) {
  if (!outcome %in% c("detect_worthy_lesions", "final_result_cat_neg")) stop("invalid outcome variable")
  ecoli_vars <- paste0("ecoli_pos", ifelse(pks_strat, "+pks_classification+","+")) 
  adj_vars <- c("senter",
                "age_invitation",
                "kjonn",
                "wcrf_index_main",
                "Smoking",
                "Utdanning",
                "colo_hemorrhoids",
                "colo_IBD",
                "antibiotics_reg_quest_comb")
  tmp_formula <- as.formula(paste0(outcome, "~", ecoli_vars, paste0(adj_vars, collapse = "+")))
  
  if (outcome == "detect_worthy_lesions") {
    tmp_mod <- glm(tmp_formula, data = dat, family = "binomial")
  } else {
    tmp_mod <- multinom(tmp_formula, data = dat)
  }
  tmp_mod %>% 
    tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
    filter(str_detect(term, "^(pks|ecoli)"))
}

ecoli_pres_abs_models <-
  tibble(outcome_var = rep(c("detect_worthy_lesions", "final_result_cat_neg"), each = 2),
         pks_strat = c(FALSE, TRUE, FALSE, TRUE)) %>% 
  mutate(cb_dat = lapply(seq(nrow(.)), function(x) cb_for_mod)) %>% 
  rowwise() %>% 
  mutate(ecoli_eval = list(evaluate_ecoli_assoc(dat = cb_dat, outcome = outcome_var, pks_strat = pks_strat))) %>% 
  select(-cb_dat) %>% 
  unnest(ecoli_eval) %>% 
  mutate(test_level = ifelse(is.na(y.level), "CRC-related findings", y.level) %>% 
           factor(levels = c("CRC-related findings", levels(meta_dat$final_result_cat_neg)[-1]))) %>% 
  mutate(pks_lab = case_when(term == "ecoli_pospresent" & pks_strat == FALSE ~ "Any E coli detected",
                             term == "ecoli_pospresent" & pks_strat == TRUE ~ "pks negative E coli detected",
                             str_detect(term, "pks positive") & pks_strat == TRUE ~ "pks positive E coli detected")) %>% 
  mutate(pks_strat = factor(pks_strat, labels = c("E coli detection overall", "E coli detection stratified by pks status")))

ecoli_pres_abs_models %>% 
  write_tsv("results/tables/ecoli_pks/ecoli_pks_detection_outcomes.tsv")