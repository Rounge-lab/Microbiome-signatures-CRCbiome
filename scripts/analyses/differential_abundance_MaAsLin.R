
differential_abundance_main_outcome <- function(models = c("main")) {
  
  tmp_meta <-
    sample_data %>% 
    select(sample_id, deltaker_id) %>% 
    left_join(meta_dat %>% select(deltaker_id, kjonn, age_invitation, senter, wcrf_index_main, 
                                  final_result, final_result_cat_neg, detect_worthy_lesions, 
                                  antibiotics_reg_quest_comb, PPI_antacids_reg_quest_comb,
                                  Utdanning, diarrhea = symtomduration_diare_cat,
                                  changing_bowel_habits = symtomduration_endrede_avforingsvaner_cat,
                                  bowel_pain = symtomduration_smerter_cat), 
              by = "deltaker_id") %>% 
    select(-c(deltaker_id))
  
  dmm <- load_dmm(synthetic_data = TRUE) %>% 
    select(sample_id, dmm = gr)
  
  identical(rownames(ko_hum), rownames(mphlan_abundance))
  cor_matrix <- cor(ko_hum[,-1], mphlan_abundance[,-1])
  
  high_cor_paths <- colnames(ko_hum[,-1])[apply(cor_matrix, 1, function(x) any(x > 0.5))]
  
  ko_hum_filtered <- 
    ko_hum %>% 
    select(-all_of(high_cor_paths))
  
  strat_dat <-
    sample_data %>% 
    left_join(dmm, by = "sample_id") %>% 
    left_join(screening_data %>% 
                rowwise() %>% 
                mutate(colo_any = ifelse(any(c_across(starts_with("colo") & !ends_with("count")) == "Yes"), "Yes", "No")) %>% 
                ungroup() %>% 
                select(deltaker_id, senter, kjonn, id, colo_diverticulitis, 
                       colo_hemorrhoids, colo_IBD, colo_any), by = "deltaker_id") %>% 
    left_join(meta_dat %>% select(deltaker_id, 
                                  antibiotics_reg_quest_comb, 
                                  PPI_antacids_reg_quest_comb, 
                                  first_round), by = "deltaker_id") %>% 
    select(sample_id, dmm, senter, kjonn, antibiotics_reg_quest_comb, 
           first_round, PPI_antacids_reg_quest_comb, 
           colo_diverticulitis, colo_hemorrhoids, colo_IBD, colo_any)
  
  
  # env ---------------------------------------------------------------------
  
  da_models <- 
    expand.grid(dataset = c("metaphlan", "keggs", "mags"),
                strat_var = c("", "dmm", "kjonn", "senter", 
                              "antibiotics_reg_quest_comb", "first_round", "PPI_antacids_reg_quest_comb",
                              "colo_diverticulitis", "colo_hemorrhoids", "colo_IBD", "colo_any"),
                model = c("unadjusted", "simplified_model", "standard_model",
                          "fully_adjusted_model", "add_dmm", "add_dna_conc", 
                          "symptoms", "antibiotics_reg_quest_comb", "PPI_antacids_reg_quest_comb"),
                outcome = c("outcome2", "outcome3"), 
                stringsAsFactors = FALSE) %>% 
    tibble() %>% 
    ## Differently adjusted models only run with unstratified data
    filter(!(model %in% c("unadjusted", "add_dmm", "add_dna_conc", "simplified_model", 
                          "fully_adjusted_model", "symptoms", 
                          "antibiotics_reg_quest_comb", "PPI_antacids_reg_quest_comb") & !strat_var %in% "")) %>% 
    bind_rows(tibble(dataset = "metaphlan",
                     strat_var = "",
                     model = "standard_model",
                     outcome = "outcome3",
                     diff_ref = levels(factor(tmp_meta$final_result_cat_neg)))) %>% 
    bind_rows(expand.grid(dataset = "metaphlan",
                          strat_var = "",
                          model = c("simplified_model",
                                    "standard_model",
                                    "fully_adjusted_model",
                                    "antibiotics_reg_quest_comb"),
                          outcome = c("lesion_loc_overall",
                                      "lesion_loc_size_distal", "lesion_loc_size_proximal",
                                      "lesion_dist_distal", "lesion_dist_proximal"),
                          diff_ref = NA,
                          stringsAsFactors = FALSE) %>% 
                tibble() %>% 
                filter(model %in% "standard_model" | outcome %in% "lesion_loc_overall")) %>% 
    bind_rows(expand.grid(dataset = c("species"),
                          strat_var = "",
                          model = c("standard_model"),
                          outcome = c("outcome2", "outcome3"),
                          diff_ref = NA,
                          stringsAsFactors = FALSE) %>% 
                tibble()) %>% 
    bind_rows(expand.grid(dataset = c("colibactin_strat"),
                          strat_var = "",
                          model = c("standard_model"),
                          outcome = c("outcome2", "outcome3"),
                          diff_ref = NA,
                          stringsAsFactors = FALSE) %>% 
                tibble()) %>% 
    mutate(res_path = paste0("data/diff_abund/Maaslin2/main_outcome/"))
  
  perform_da <- function(dataset, strat_var, model, outcome, diff_ref = NA, prev_level = 0.1, res_path) {
    
    if (dataset == "metaphlan") tmp_abund <- mphlan_abundance
    if (dataset == "keggs") tmp_abund <- ko_hum_filtered
    if (dataset == "mags") tmp_abund <- MAGs_abundance
    if (dataset == "colibactin_strat") {
      tmp_abund <-
        mphlan_abundance %>% 
        left_join(read_tsv("data/input_processed/pks_detection.tsv", col_types = cols()) %>% 
                    select(sample_id, pks_classification), by = "sample_id") %>% 
        mutate(pks_pos_ecoli = case_when(str_detect(pks_classification, "positive") ~ SGB10068,
                                         TRUE ~ 0),
               pks_neg_ecoli = case_when(str_detect(pks_classification, "negative") ~ SGB10068,
                                         TRUE ~ 0)) %>% 
        select(-SGB10068)
    }  
    
    res_path <- paste0(res_path, "/", dataset, "/", outcome, "/", model)
    
    if (strat_var %in% "" & !is.na(diff_ref)) {
        res_path <- paste0(res_path, "/pairwise")
    }
    if (!strat_var %in% "") {
      res_path <- paste0(res_path, "/strat/", strat_var)
    }
    if (prev_level != 0.1) {
      res_path <- paste0(res_path, "/low_prev")
    }
    
    if (!dir.exists(res_path)) dir.create(res_path, recursive = TRUE)
    
    if (model %in% c("fully_adjusted_model", "add_dna_conc", "add_dmm", "symptoms", "antibiotics_reg_quest_comb", "PPI_antacids_reg_quest_comb")) {
      metadata <- 
        sample_data %>% 
        select(sample_id, deltaker_id, Total_Reads_QC_ATLAS) %>% 
        left_join(screening_data %>% 
                    select(deltaker_id, id, colo_hemorrhoids), by = "deltaker_id") %>% 
        left_join(meta_dat %>% 
                    select(deltaker_id, senter, age_invitation, kjonn, PhysAct_Score, BMI, Alko, 
                           Arbeid_lump, Bowel_disorder_merged,
                           Smoking, Utdanning, final_result, detect_worthy_lesions, final_result_cat_neg), by = "deltaker_id") %>% 
        select(-c(id, deltaker_id))
      if (model %in% "add_dna_conc") {
        metadata <- metadata %>% left_join(sample_data %>% select(sample_id, qubit_total_dna), by = "sample_id")
      }
      if (model %in% "add_dmm") {
        metadata <- metadata %>% left_join(dmm, by = "sample_id")
      }
      if (model %in% "symptoms") {
        metadata <- metadata %>% left_join(tmp_meta %>% select(sample_id, diarrhea, changing_bowel_habits, bowel_pain), by = "sample_id")
      }
      if (model %in% "antibiotics_reg_quest_comb") {
        metadata <- metadata %>% left_join(strat_dat %>% select(sample_id, antibiotics_reg_quest_comb), by = "sample_id")
      }
      if (model %in% "PPI_antacids_reg_quest_comb") {
        metadata <- metadata %>% left_join(tmp_meta %>% select(sample_id, PPI_antacids_reg_quest_comb), by = "sample_id") 
      }
    } else if (model %in% c("simplified_model", "standard_model")) {
      metadata <- 
        sample_data %>% 
        select(sample_id, deltaker_id, Total_Reads_QC_ATLAS) %>% 
        left_join(meta_dat %>% 
                    select(deltaker_id, senter, age_invitation, kjonn, wcrf_index_main, Smoking, Utdanning,
                           final_result, detect_worthy_lesions, final_result_cat_neg), by = "deltaker_id") %>% 
        select(-c(deltaker_id))
      if (model %in% "standard_model") {
        metadata <- 
          metadata %>% 
          left_join(sample_data %>% select(sample_id, deltaker_id), by = join_by(sample_id)) %>% 
          left_join(screening_data %>% select(deltaker_id, colo_hemorrhoids, colo_IBD), by = join_by(deltaker_id)) %>% 
          left_join(tmp_meta %>% select(sample_id, antibiotics_reg_quest_comb), by = join_by(sample_id)) %>% 
          # mutate(across(.cols = starts_with("colo"), .fn = function(x) x %>% fct_recode(a_No = "No"))) %>% 
          select(-deltaker_id)
      }
    }
    if (model %in% "unadjusted") {
      metadata <- 
        sample_data %>% 
        select(sample_id, deltaker_id) %>% 
        left_join(meta_dat %>% 
                    select(deltaker_id, final_result, final_result_cat_neg, detect_worthy_lesions), by = "deltaker_id") %>% 
        select(-c(deltaker_id))
    }
    
    adj_mod <- 
      metadata %>% 
      (function(x) {
        if (strat_var == "colo_any") {
          x %>% 
            select(-starts_with("colo_"))
        } else {
          x
        }
      }) %>% 
      select(-c(final_result, detect_worthy_lesions, final_result_cat_neg, sample_id, any_of(strat_var))) %>% 
      names() %>% 
      paste(collapse = ",")
    
    adj_mod <- ifelse(nchar(adj_mod) > 0, paste0(",",adj_mod), "")
    
    variables <- bind_rows(variables, tibble(var_id = "dmm", ref = "V1")) %>% distinct()
    ref <- metadata %>% select(-any_of(strat_var)) %>% names() %>% enframe(value = "var_id") %>% inner_join(variables %>% select(var_id, ref), by = "var_id") %>% filter(!is.na(ref))
    
    if (strat_var == "") {
      metadata <-
        metadata %>% 
        mutate(strat = "no_strat")
    } else if (!strat_var %in% names(metadata)) {
      metadata <-
        metadata %>% 
        left_join(strat_dat %>% select(sample_id, any_of(strat_var)), by = "sample_id") %>% 
        select(everything(), strat = all_of(strat_var))
    } else {
      metadata <-
        metadata %>% 
        select(everything(), strat = all_of(strat_var))
    }
    
    if (outcome == "outcome1") {
      tmp_mod <- paste0("final_result", adj_mod)
      if (is.na(diff_ref)) {
        ref <- ref  %>% 
          # filter(var_id != "detect_worthy_lesions") %>% 
          filter(!str_detect(var_id, "(final_result_cat_neg|detect_worthy_lesions)")) %>% 
          mutate(ref = paste(var_id, ref, sep = ",")) %>% pull(ref) %>% paste(collapse = ";")  
      } else {
        ref <- ref  %>% 
          # filter(var_id != "detect_worthy_lesions") %>% 
          filter(!str_detect(var_id, "(final_result_cat_neg|detect_worthy_lesions)")) %>% 
          mutate(ref = case_when(var_id %in% "final_result" ~ paste0("0", diff_ref),
                                 TRUE ~ ref)) %>% 
          mutate(ref = paste(var_id, ref, sep = ",")) %>% pull(ref) %>% paste(collapse = ";")  
      }
    } else if (outcome == "outcome2") {
      tmp_mod <- paste0("detect_worthy_lesions", adj_mod)
      ref <- ref  %>% 
        filter(!str_detect(var_id, "^final_result")) %>% 
        mutate(ref = paste(var_id, ref, sep = ",")) %>% pull(ref) %>% paste(collapse = ";")
    } else if (outcome == "outcome3") {
      tmp_mod <- paste0("final_result_cat_neg", adj_mod)
      if (is.na(diff_ref)) {
        ref <- ref  %>% 
          # filter(var_id != "detect_worthy_lesions") %>% 
          filter(!str_detect(var_id, "(final_result$|detect_worthy_lesions)")) %>% 
          mutate(ref = paste(var_id, ref, sep = ",")) %>% pull(ref) %>% paste(collapse = ";")  
      } else {
        ref <- ref  %>% 
          # filter(var_id != "detect_worthy_lesions") %>% 
          filter(!str_detect(var_id, "(final_result$|detect_worthy_lesions)")) %>% 
          mutate(ref = case_when(var_id %in% "final_result_cat_neg" ~ paste0("0", diff_ref),
                                 TRUE ~ ref)) %>% 
          mutate(ref = paste(var_id, ref, sep = ",")) %>% pull(ref) %>% paste(collapse = ";")  
      }
    }
    
    if (str_detect(outcome, "^lesion_")) {
      if (outcome %in% "lesion_loc_overall") {
        location_vars <- c("distal_acn", "proximal_acn")
      } else if (outcome %in% "lesion_loc_size_distal") {
        location_vars <- "sum_lesion_diameter_distal"
      } else if (outcome %in% "lesion_loc_size_proximal") {
        location_vars <- c("sum_lesion_diameter_proximal")
      } else if (outcome %in% "lesion_dist_distal") {
        location_vars <- c("most_serious_distancefromanus_processed")
        filter_side <- "Distal"
      } else if (outcome %in% "lesion_dist_proximal") {
        location_vars <- c("most_serious_distancefromanus_processed")
        filter_side <- "Proximal"
      }
      
      tmp_mod <- paste0(paste0(location_vars, collapse = ","), adj_mod)
      
      ref <- ref  %>% 
        filter(!var_id %in% c("final_result", "detect_worthy_lesions", "final_result_cat_neg")) %>% 
        mutate(ref = paste(var_id, ref, sep = ",")) %>% pull(ref) %>% paste(collapse = ";")
      
      metadata <- 
        metadata %>% 
        left_join(sample_data %>% select(sample_id, deltaker_id), by = join_by(sample_id)) %>% 
        left_join(screening_data %>% select(deltaker_id, all_of(location_vars)), by = join_by(deltaker_id))
      
      if (outcome %in% c("lesion_loc_size_distal", "lesion_loc_size_distal")) {
        metadata <- 
          metadata %>% 
          filter(!!sym(location_vars[1]) > 0)
      }
      
      if (outcome %in% c("lesion_dist_distal", "lesion_dist_proximal")) {
        metadata <- 
          metadata %>% 
          inner_join(screening_data %>% 
                      inner_join(lesion_locs %>% 
                                   select(most_serious_segment = segment, side_processed) %>% 
                                   distinct() %>%
                                   filter(side_processed == filter_side) %>% 
                                   na.omit(), by = join_by(most_serious_segment)) %>% 
                      select(deltaker_id),
                    by = join_by(deltaker_id))
      }
    }
    
    metadata %>% 
      group_by(strat) %>% 
      group_split() %>% 
      lapply(function(metadata_strat) {
        
        ## Add folder for stratification
        if (!metadata_strat$strat[1] == "no_strat") {
          res_path <- paste0(res_path, "/", metadata_strat$strat[1])
          if (!dir.exists(res_path)) dir.create(res_path, recursive = TRUE)
        }
        ## Add folder for reference level if different reference level is set
        if (!is.na(diff_ref)) {
          shorthand_fin_res <- tibble(outcome = levels(factor(metadata$final_result_cat_neg)),
                                      sh = c("neg", "naa", "as", "aa", "crc"))
          res_path <- paste0(res_path, "/ref_", shorthand_fin_res %>% filter(outcome %in% diff_ref) %>% pull(sh))
          if (!dir.exists(res_path)) dir.create(res_path, recursive = TRUE)
        }
        
        tmp_maaslin <-
          Maaslin2::Maaslin2(
            input_data = metadata_strat %>% 
              select(sample_id) %>% 
              left_join(tmp_abund) %>% 
              as.data.frame() %>% 
              column_to_rownames("sample_id"),
            input_metadata = metadata_strat %>% 
              (function(x) if (!is.na(diff_ref)) x %>% mutate(final_result_cat_neg = ifelse(final_result_cat_neg %in% diff_ref, 
                                                                                            paste0("0", as.character(final_result_cat_neg)), 
                                                                                            paste0("1", as.character(final_result_cat_neg)))) else x) %>% 
              as.data.frame() %>% 
              column_to_rownames("sample_id"),
            min_prevalence = prev_level,
            transform = "LOG",
            normalization = "TSS",
            analysis_method = "LM",
            output = res_path,
            fixed_effects = tmp_mod,
            reference = ref,
            plot_heatmap = FALSE,
            plot_scatter = FALSE
          )
        
        tmp_maaslin$results %>% 
          tibble() %>% 
          mutate(model = tmp_mod,
                 dataset = dataset,
                 strat = metadata_strat$strat[1]) %>% 
          (function(x) if (!is.na(diff_ref)) x %>% mutate(ref = diff_ref) else x) %>% 
          (function(x) {
            if (str_detect(outcome, "^lesion_")) {
              x %>% filter(metadata %in% location_vars)
            } else {
              x %>% filter(str_detect(metadata, "(final_result|detect_worthy_lesions|final_result_cat_neg)"))
            }
          })
      }) %>% 
      bind_rows()
  }
  
  if ("main" %in% models) {
    ## Standard models
    diff_abund_main_results <-
      da_models %>%
      filter(strat_var %in% "",
             model %in% "standard_model",
             dataset %in% c("metaphlan","keggs", "mags"),
             is.na(diff_ref),
             outcome %in% c("outcome2", "outcome3")) %>% 
      rowwise() %>% 
      mutate(da_res = list(perform_da(dataset, strat_var, model, outcome, 
                                      res_path = "data/diff_abund/Maaslin2/main_outcome"))) %>% 
      ungroup()
    
    diff_abund_main_results %>% 
      select(da_res) %>%
      unnest(da_res) %>% 
      write_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_main_models.tsv")  
  }
  
  if ("strat" %in% models) {
    ## Stratified models
    diff_abund_strat_results <-
      da_models %>% 
      filter(!strat_var %in% "") %>% 
      rowwise() %>%
      mutate(da_res = list(perform_da(dataset, strat_var, model, outcome, 
                                      res_path = "data/diff_abund/Maaslin2/main_outcome"))) %>% 
      ungroup()
    
    diff_abund_strat_results %>% 
      select(da_res, strat_var) %>%
      unnest(da_res)%>% 
      write_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_strat_models.tsv")  
  }
  
  if ("adj" %in% models) {
    ## Adjusted models
    diff_abund_adj_res <-
      da_models %>% 
      filter(strat_var %in% "",
             !model %in% c("unadjusted", "simplified_model", "standard_model"),
             outcome %in% c("outcome2", "outcome3")) %>% 
      rowwise() %>% 
      mutate(da_res = list(perform_da(dataset, strat_var, model, outcome, 
                                      res_path = "data/diff_abund/Maaslin2/main_outcome"))) %>% 
      ungroup()
    
    diff_abund_adj_res %>% 
      select(da_res) %>%
      unnest(da_res) %>% 
      write_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_full_adj_models.tsv")
  }
  if ("pairwise" %in% models) {
    diff_abund_pairwise_da_models <-
      da_models %>% 
      filter(!is.na(diff_ref)) %>% 
      rowwise() %>% 
      mutate(da_res = list(perform_da(dataset, strat_var, model, outcome, diff_ref, 
                                      res_path = "data/diff_abund/Maaslin2/main_outcome"))) %>% 
      ungroup()
    
    diff_abund_pairwise_da_models %>% 
      select(da_res) %>%
      unnest(da_res) %>% 
      write_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_pairwise_models.tsv")
  }
  if ("low_prev" %in% models) {
    diff_abund_low_prev_da_models <-
      da_models %>% 
      filter(is.na(diff_ref)) %>%
      filter(strat_var %in% "") %>% 
      filter(model %in% c("standard_model")) %>% 
      filter(dataset %in% c("metaphlan", "species")) %>% 
      filter(str_detect(outcome, "outcome")) %>% 
      rowwise() %>% 
      mutate(da_res = list(perform_da(dataset, strat_var, model, outcome, diff_ref, prev_level = 0.001, 
                                      res_path = "data/diff_abund/Maaslin2/main_outcome"))) %>% 
      ungroup()
    
    diff_abund_low_prev_da_models %>% 
      select(da_res) %>%
      unnest(da_res) %>% 
      write_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_low_prev_models.tsv")
  }
  if ("unadjusted" %in% models) {
    diff_abund_unadjusted_da_models <-
      da_models %>% 
      filter(model %in% c("unadjusted")) %>% 
      rowwise() %>% 
      mutate(da_res = list(perform_da(dataset, strat_var, model, outcome, diff_ref, prev_level = 0.1, 
                                      res_path = "data/diff_abund/Maaslin2/main_outcome"))) %>% 
      ungroup()
    
    diff_abund_unadjusted_da_models %>% 
      select(da_res) %>%
      unnest(da_res) %>% 
      write_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_unadjusted_models.tsv")
  }
  if ("loc" %in% models) {
    diff_abund_loc_da_models <-
      da_models %>% 
      filter(str_detect(outcome, "^lesion")) %>% 
      rowwise() %>% 
      mutate(da_res = list(perform_da(dataset, strat_var, model, 
                                      outcome, diff_ref, prev_level = 0.1, 
                                      res_path = "data/diff_abund/Maaslin2/main_outcome"))) %>% 
      ungroup()
    
    diff_abund_loc_da_models %>% 
      select(da_res, outcome) %>%
      unnest(da_res) %>% 
      write_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_localization_models.tsv")
  }
  if ("colibactin_strat" %in% models) {
    diff_abund_colibactin_strat_da_models <-
      da_models %>% 
      filter(str_detect(dataset, "colibactin_strat")) %>% 
      rowwise() %>% 
      mutate(da_res = list(perform_da(dataset, strat_var, model, 
                                      outcome, diff_ref, prev_level = 0.1, 
                                      res_path = "data/diff_abund/Maaslin2/main_outcome"))) %>% 
      ungroup()
    
    diff_abund_colibactin_strat_da_models %>% 
      select(da_res, outcome) %>%
      unnest(da_res) %>% 
      write_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_colibactin_strat_models.tsv")
  }
}

