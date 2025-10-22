

beta_by_level <- function(dists, m_dat, test_var, adjust_vars, ref_level_test_var, id_col = "sample_id") {
  
  tmp_meta <-
    dists %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column("id") %>% 
    tibble() %>% 
    select(id) %>% 
    inner_join(m_dat %>% 
                 select(all_of(unique(c(id_col, test_var, adjust_vars)))) %>% 
                 rename(id = all_of(id_col)),
               by = "id")
  
  if (!tmp_meta %>% 
      pull(all_of(test_var)) %>% 
      is.double()) {
    
    tmp_lvls <-
      tmp_meta %>% 
      pull(all_of(test_var)) %>% 
      factor() %>% 
      relevel(ref_level_test_var)
    
    lapply(levels(tmp_lvls)[-1], function(x) {
      tmp_sel <- tmp_lvls %in% c(ref_level_test_var, x) &
        !apply(is.na(tmp_meta[,-1]), 1, any) 
      
      tmp <- vegan::adonis2(formula = as.dist(as.matrix(dists)[tmp_sel, tmp_sel]) ~ ., data = tmp_meta %>% filter(tmp_sel) %>% select(-id), by = "margin", permutations = ifelse(synthetic_data, 10, 999))
      tmp2 <- adonis_OmegaSq(tmp) %>% as.data.frame() %>% rownames_to_column("term") %>% tibble() %>% bind_cols("R2" = tmp$R2)
      
      tmp2 %>% 
        mutate(test_var = test_var,
               ref_level = ref_level_test_var,
               test_level = x,
               n = sum(tmp_sel)) %>% 
        filter(!term %in% c("Residual","Total"))
    }) %>% bind_rows()
  } else {
    tmp_sel <- !apply(is.na(tmp_meta[,-1]), 1, any) 
    
    tmp <- vegan::adonis2(formula = as.dist(as.matrix(dists)[tmp_sel, tmp_sel]) ~ ., data = tmp_meta %>% filter(tmp_sel) %>% select(-id), by = "margin", permutations = ifelse(synthetic_data, 10, 999))
    tmp2 <- adonis_OmegaSq(tmp) %>% as.data.frame() %>% rownames_to_column("term") %>% tibble() %>% bind_cols("R2" = tmp$R2)
    
    tmp2 %>% 
      mutate(test_var = test_var,
             ref_level = NA,
             test_level = NA,
             n = sum(tmp_sel)) %>% 
      filter(!term %in% c("Residual","Total"))
  }
  
}

get_strings_from_vec <- function(strings, vec) {
  sapply(strings, function(x) {
    str_extract(vec, x)
  }) %>% 
    apply(1, function(x) {
      if (any(!is.na(x))) {
        x[which(!is.na(x))]  
      } else {
        NA
      }
    }) 
}

tmp_dats <- c("mags", "metaphlan", "ko_hum")

alpha_diversity <- 
  lapply(tmp_dats, function(i_dat) {
    dat <- list("mags" = MAGs_abundance, 
                "metaphlan" = mphlan_abundance, 
                "ko_hum" = ko_hum)[[i_dat]]
    
    print(paste0("dataset ", i_dat, ": ", which(tmp_dats %in% i_dat), " of ", length(tmp_dats)))
    
    dat %>% 
      pivot_longer(-sample_id) %>% 
      group_by(sample_id) %>% 
      summarize(observed = sum(value > 0),
                shannon = diversity(value),
                invsimpson = diversity(value, index = "invsimpson")) %>% 
      ungroup() %>% 
      mutate(dataset = i_dat)
  }) %>% 
  bind_rows()

alpha_cat_tests <-
  lapply(as.character(variables$var_id), function(i) {
    
    print(paste0("variable ", i, ": ", which(variables$var_id %in% i), " of ", length(variables$var_id)))
    
    ref_lvl <- variables %>% filter(var_id %in% i) %>% pull(ref)
    
    alpha_diversity %>%
      pivot_longer(c(observed, shannon, invsimpson), names_to = "index", values_to = "diversity") %>% 
      group_by(dataset, index) %>% 
      group_split() %>% 
      lapply(function(x) {
        tmp <-
          x %>% 
          left_join(meta_dat_cat %>% 
                      left_join(sample_data %>% 
                                  select(sample_id, deltaker_id, Total_Reads_QC_ATLAS),
                                by = "deltaker_id") %>% 
                      select(sample_id, any_of(c("kjonn", "age_invitation", "senter", "Total_Reads_QC_ATLAS", i))), 
                    by = "sample_id") %>% 
          select(-c(dataset, index))
        
        if (is.na(ref_lvl)) {
          ref_lvl <- "low"
        }
        
        if (!is.na(ref_lvl)) {
          tmp <-
            tmp %>% 
            mutate(across(.cols = matches(i), .fns = function(s) factor(s) %>% relevel(ref = ref_lvl)))
        }
        
        tmp %>% 
          select(-sample_id) %>% 
          lm(formula = diversity ~ ., data = .) %>% 
          tidy(conf.int = TRUE) %>%
          mutate(variable = get_strings_from_vec(names(tmp)[-c(1,2)], term),
                 lvl = str_remove(term, variable),
                 index = x$index[1],
                 dataset = x$dataset[1],
                 test_var = i,
                 ref_level = ref_lvl)
      }) %>% 
      bind_rows()
  }) %>% 
  bind_rows()

# if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
# 
# alpha_cat_tests %>% 
#   write_tsv("results/tables/alpha_host_vars_cat.tsv"))

permanova_cat <- 
  lapply(tmp_dats, function(i_dat) {
    dat <- list("mags" = MAGs_abundance, 
                "metaphlan" = mphlan_abundance, 
                "ko_hum" = ko_hum)[[i_dat]]
    
    print(paste0("dataset ", i_dat, ": ", which(tmp_dats %in% i_dat), " of ", length(tmp_dats)))
    
    tmp_dists <- 
      dat %>% 
      as.data.frame() %>% 
      column_to_rownames("sample_id") %>% 
      as.matrix() %>% 
      vegdist(method = "bray")
    
    lapply(as.character(variables$var_id), function(i) {
      
      print(paste0("variable ", i, ": ", which(variables$var_id %in% i), " of ", length(variables$var_id)))
      
      ref_lvl <- variables %>% filter(var_id %in% i) %>% pull(ref)
      
      tmp_dists %>%
        beta_by_level(m_dat = meta_dat_cat %>% left_join(sample_data %>% select(sample_id, deltaker_id, Total_Reads_QC_ATLAS), by = "deltaker_id"),
                      test_var = i,
                      adjust_vars = c("kjonn", "age_invitation", "senter", "Total_Reads_QC_ATLAS"),
                      ref_level_test_var = ifelse(is.na(ref_lvl), "low", ref_lvl))
    }) %>% 
      bind_rows() %>% 
      mutate(dataset = i_dat)
  }) %>% 
  bind_rows() %>% 
  tibble()

# if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
# 
# permanova_cat %>%
#   write_tsv("results/tables/permanova_host_vars_cat.tsv")


# multinom models ---------------------------------------------------------

alpha_adj_outcome_tests_multinom <-
  lapply(c("standard_model"), function(mod) {
    
    lapply(as.character(c("detect_worthy_lesions", "final_result_cat_neg")), function(i) {
      
      if (mod %in% "standard_model") {
        adj_vars <- 
          sample_data %>% 
          select(sample_id, deltaker_id, Total_Reads_QC_ATLAS) %>% 
          left_join(meta_dat %>% 
                      select(deltaker_id, senter, age_invitation, kjonn, wcrf_index_main, 
                             Smoking, Utdanning, colo_hemorrhoids, colo_IBD, 
                             antibiotics_reg_quest_comb, any_of(i)), by = "deltaker_id") %>% 
          select(-c(deltaker_id))
      }
      
      lapply(c("all", 
               "no_recent_antibiotics", 
               "Male", "Female", 
               "Moss", "Bærum",
               "proximal_acn", "distal_acn"), function(excl_cat) {
                 
                 if (excl_cat == "no_recent_antibiotics") {
                   adj_vars <- 
                     adj_vars %>% 
                     inner_join(meta_dat %>% 
                                  filter(!antibiotics_reg_quest_comb == "Yes") %>% 
                                  left_join(sample_data %>% select(sample_id, deltaker_id), by = "deltaker_id") %>% 
                                  select(sample_id), by = "sample_id") %>% 
                     select(-any_of("antibiotics_reg_quest_comb"))
                 }
                 if (excl_cat %in% c("Male", "Female")) {
                   adj_vars <- adj_vars %>% filter(kjonn %in% excl_cat) %>% select(-kjonn)
                 }
                 if (excl_cat %in% c("Bærum", "Moss")) {
                   adj_vars <- adj_vars %>% filter(senter %in% excl_cat) %>% select(-senter)
                 }
                 if (excl_cat %in% c("proximal_acn", "distal_acn")) {
                   adj_vars %>% 
                     inner_join(screening_data %>% 
                                  rename(tmp = any_of("distal_acn")) %>% 
                                  filter(tmp %in% 1 | detect_worthy_lesions %in% "negative") %>% 
                                  left_join(sample_data, by = "deltaker_id") %>% 
                                  select(sample_id), by = "sample_id")
                 }
                 
                 print(paste0("variable ", i, ": ", which(variables$var_id %in% i), " of 2; subset: ", excl_cat))
                 
                 ref_lvl <- variables %>% filter(var_id %in% i) %>% pull(ref)
                 
                 alpha_diversity %>%
                   inner_join(adj_vars, by = "sample_id") %>% 
                   rename(dep_var = any_of(i)) %>% 
                   pivot_longer(c(observed, shannon, invsimpson), names_to = "index", values_to = "diversity") %>% 
                   group_by(dataset, index) %>% 
                   group_split() %>% 
                   lapply(function(x) {
                     tmp <-
                       x %>% 
                       select(-c(dataset, index))
                     
                     if (is.na(ref_lvl)) {
                       ref_lvl <- "low"
                     }
                     
                     if (!is.na(ref_lvl)) {
                       tmp <-
                         tmp %>% 
                         mutate(across(.cols = dep_var, .fns = function(s) factor(s) %>% relevel(ref = ref_lvl)))
                     } 
                     
                     tmp %>% 
                       select(-sample_id) %>% 
                       mutate(across(where(is.numeric), ~ c(scale(.)))) %>% 
                       multinom(dep_var ~ ., data = .) %>% 
                       tidy(conf.int = TRUE, exponentiate = TRUE) %>%
                       mutate(variable = get_strings_from_vec(names(tmp)[-c(1,2)], term),
                              lvl = str_remove(term, variable),
                              index = x$index[1],
                              dataset = x$dataset[1],
                              test_var = i,
                              ref_level = ref_lvl,
                              subset_cat = excl_cat,
                              model = mod)
                   }) %>% 
                   bind_rows()
               }) %>% 
        bind_rows()
    }) %>% 
      bind_rows()
  }) %>% 
  bind_rows() %>% 
  ## Proximal/distal should not be evaluated for final_result
  filter(!(test_var %in% "final_result_cat_neg" & str_detect(subset_cat, "acn")))

# if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
# 
# alpha_adj_outcome_tests_multinom %>% 
#   write_tsv("results/tables/alpha_adj_outcome_tests.tsv"))


dist_dat <-
  mphlan_abundance %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "sample_id") %>% 
  as.matrix() %>% 
  vegan::vegdist()


fin_res_permanova <- 
  lapply(c("detect_worthy_lesions", "final_result_cat_neg"), function(i) {
    
    lapply(c("unadjusted", "standard_model"), function(mod) {
      
      if (mod %in% "unadjusted") {
        adj_vars <- 
          sample_data %>% 
          select(sample_id, deltaker_id) %>% 
          left_join(meta_dat %>% 
                      select(deltaker_id, any_of(i)), by = "deltaker_id") %>% 
          select(-c(deltaker_id))
      }
      if (mod %in% "standard_model") {
        adj_vars <- 
          sample_data %>% 
          select(sample_id, deltaker_id, Total_Reads_QC_ATLAS) %>% 
          left_join(meta_dat %>% 
                      select(deltaker_id, senter, age_invitation, kjonn, wcrf_index_main, 
                             Smoking, Utdanning, colo_hemorrhoids, colo_IBD, 
                             antibiotics_reg_quest_comb, any_of(i)), by = "deltaker_id") %>% 
          select(-c(deltaker_id))
      }
      
      tmp_meta <- 
        adj_vars %>% 
        na.omit()
      
      tmp_sel <- colnames(as.matrix(dist_dat)) %in% tmp_meta$sample_id
      
      tmp_meta <- 
        tmp_meta %>% 
        right_join(colnames(as.matrix(dist_dat))[tmp_sel] %>% 
                     enframe(value = "sample_id") %>% 
                     select(sample_id), by = join_by(sample_id)) %>% 
        select(-sample_id)
      set.seed(1)
      tmp <- vegan::adonis2(formula = as.dist(as.matrix(dist_dat)[tmp_sel, tmp_sel]) ~ ., data = tmp_meta, by = "margin")
      tmp2 <- adonis_OmegaSq(tmp) %>% as.data.frame() %>% rownames_to_column("term") %>% tibble() %>% bind_cols("R2" = tmp$R2)
      
      tmp2 %>% 
        mutate(adjustment = mod) %>% 
        filter(!term %in% c("Residual","Total"))
    }) %>% 
      bind_rows()
  }) %>% 
  bind_rows()

# if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)

# fin_res_permanova %>% 
#   write_tsv("results/tables/beta_outcome_tests.tsv")


# Test pairwise differences between outcome groups ------------------------

## Pairwise differences will be tested only in the standard model

alpha_adj_pairwise_outcome_tests_multinom <-
  lapply(unique(meta_dat$final_result_cat_neg), function(ref_outcome) {
      
      mod <- "standard_model"
      i <- "final_result_cat_neg"

      if (mod %in% "standard_model") {
        adj_vars <- 
          sample_data %>% 
          select(sample_id, deltaker_id, Total_Reads_QC_ATLAS) %>% 
          left_join(meta_dat %>% 
                      select(deltaker_id, senter, age_invitation, kjonn, wcrf_index_main, 
                             Smoking, Utdanning, colo_hemorrhoids, colo_IBD, 
                             antibiotics_reg_quest_comb, any_of(i)), by = "deltaker_id") %>% 
          select(-c(deltaker_id))
      }
      
      print(paste0("Testing ", ref_outcome, "."))
      
      ref_lvl <- ref_outcome
      
      alpha_diversity %>%
        inner_join(adj_vars, by = "sample_id") %>% 
        rename(dep_var = any_of(i)) %>% 
        pivot_longer(c(observed, shannon, invsimpson), names_to = "index", values_to = "diversity") %>% 
        group_by(dataset, index) %>% 
        group_split() %>% 
        lapply(function(x) {
          tmp <-
            x %>% 
            select(-c(dataset, index)) %>% 
            mutate(across(.cols = dep_var, .fns = function(s) s %>% fct_relevel(as.character(ref_lvl))))# factor(s) %>% relevel(ref = ref_lvl)))
          
          tmp %>% 
            select(-sample_id) %>% 
            mutate(across(where(is.numeric), ~ c(scale(.)))) %>% 
            multinom(dep_var ~ ., data = .) %>% 
            tidy(conf.int = TRUE, exponentiate = TRUE) %>%
            mutate(variable = get_strings_from_vec(names(tmp)[-c(1,2)], term),
                   lvl = str_remove(term, variable),
                   index = x$index[1],
                   dataset = x$dataset[1],
                   test_var = i,
                   ref_level = ref_lvl,
                   model = mod)
        }) %>% 
        bind_rows()
  }) %>% 
  bind_rows()

# alpha_adj_pairwise_outcome_tests_multinom %>% 
#   write_tsv("results/tables/alpha_adj_pairwise_fin_res_outcome_tests.tsv")


beta_pairwise_fin_res <- function() {
  
  dist_dat <-
    mphlan_abundance %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "sample_id") %>% 
    as.matrix() %>% 
    vegan::vegdist()
  
  i <- "final_result_cat_neg"

  ## Perform pairwise permanova tests to determine group-wise differences
  pairwise_fin_res_dist_comp <-
    lapply(c("unadjusted", "standard_model"), function(mod) {
      
      if (mod %in% "unadjusted") {
        adj_vars <- 
          sample_data %>% 
          select(sample_id, deltaker_id) %>% 
          left_join(meta_dat %>% 
                      select(deltaker_id, any_of(i)), by = "deltaker_id") %>% 
          select(-c(deltaker_id))
      } else if (mod %in% "standard_model") {
        adj_vars <- 
          sample_data %>% 
          select(sample_id, deltaker_id, Total_Reads_QC_ATLAS) %>% 
          left_join(meta_dat %>% 
                      select(deltaker_id, senter, age_invitation, kjonn, wcrf_index_main, 
                             Smoking, Utdanning, colo_hemorrhoids, colo_IBD, 
                             antibiotics_reg_quest_comb, any_of(i)), by = "deltaker_id") %>% 
          select(-c(deltaker_id))
        }
      
      lapply(unique(meta_dat_cat$final_result_cat_neg), function(ref_lvl) {
        lapply(unique(meta_dat_cat$final_result_cat_neg), function(test_lvl) {
          if (ref_lvl != test_lvl) {
            
            print(paste("Adjustment:", mod, "\nReference level:", ref_lvl, "\nTest level:", test_lvl))
            
            tmp_meta <- 
              adj_vars %>% 
              filter(final_result_cat_neg %in% c(ref_lvl, test_lvl)) %>% 
              na.omit()
            
            tmp_sel <- colnames(as.matrix(dist_dat)) %in% tmp_meta$sample_id
            
            tmp_meta <- 
              tmp_meta %>% 
              right_join(colnames(as.matrix(dist_dat))[tmp_sel] %>% 
                           enframe(value = "sample_id") %>% 
                           select(sample_id), by = join_by(sample_id)) %>% 
              select(-sample_id)
            
            tmp <- vegan::adonis2(formula = as.dist(as.matrix(dist_dat)[tmp_sel, tmp_sel]) ~ ., data = tmp_meta, by = "margin")
            tmp2 <- adonis_OmegaSq(tmp) %>% as.data.frame() %>% rownames_to_column("term") %>% tibble() %>% bind_cols("R2" = tmp$R2)
            
            tmp2 %>%
              mutate(ref_level = ref_lvl,
                     test_level = test_lvl,
                     n_ref = sum(tmp_meta$final_result_cat_neg %in% ref_lvl),
                     n_test = sum(tmp_meta$final_result_cat_neg %in% test_lvl),
                     adjustment = mod) %>%
              filter(!term %in% c("Residual","Total"))
          }
        }) %>% 
          bind_rows()
      }) %>% 
        bind_rows()
    }) %>% 
    bind_rows()
  
  pairwise_fin_res_dist_comp 
}

beta_pairwise_fin_res() #%>% 
  # write_tsv("results/tables/beta_fin_res_pairwise.tsv")



# Site specificity --------------------------------------------------------



test_alpha_outcome_assoc_by_site_spec <- function(site_coding = "lesion_type_loc_dummy") {
  
  adj_vars <- 
    sample_data %>% 
    select(sample_id, deltaker_id, Total_Reads_QC_ATLAS) %>% 
    left_join(meta_dat %>% 
                select(deltaker_id, senter, age_invitation, 
                       kjonn, wcrf_index_main, Smoking, Utdanning,
                       colo_hemorrhoids, colo_IBD, antibiotics_reg_quest_comb), 
              by = "deltaker_id")
  
  if (site_coding == "lesion_type_loc_dummy") {
    adj_vars <- 
      adj_vars %>% 
      left_join(screening_data %>% 
                  select(contains("distal"), contains("proximal"), 
                         -contains("distance"), -contains("acn"), -contains("ssp_ssa"),
                         -starts_with("sum"),
                         deltaker_id), by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_loc_overall") {
    adj_vars <-
      adj_vars %>% 
      left_join(screening_data %>% 
                  select(contains("distal") & contains("acn"),
                         contains("proximal") & contains("acn"),
                         deltaker_id), by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_loc_size") {
    adj_vars <-
      adj_vars %>% 
      inner_join(screening_data %>% 
                  filter(!str_detect(final_result, "(1.|3.|4.)")) %>% 
                  select(contains("distal") & starts_with("sum"),
                         contains("proximal") & starts_with("sum"),
                         -contains("volume"), deltaker_id), 
                by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_loc_size_distal") {
    adj_vars <-
      adj_vars %>% 
      inner_join(screening_data %>%
                  filter(sum_lesion_diameter_distal > 0) %>%
                  select(sum_lesion_diameter_distal,
                         deltaker_id),
                by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_loc_size_proximal") {
    adj_vars <-
      adj_vars %>% 
      inner_join(screening_data %>%
                  filter(sum_lesion_diameter_proximal > 0) %>%
                  select(sum_lesion_diameter_proximal,
                         deltaker_id),
                by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_dist_distal") {
    adj_vars <- 
      adj_vars %>% 
      left_join(screening_data %>% 
                  inner_join(lesion_locs %>% 
                               select(most_serious_segment = segment, side_processed) %>% 
                               distinct() %>%
                               filter(side_processed == "Distal") %>% 
                               na.omit(), by = join_by(most_serious_segment)) %>% 
                  select(most_serious_distancefromanus_processed, deltaker_id, final_result),
                by = join_by(deltaker_id)) %>% 
      select(-c(deltaker_id, final_result))
  } else if (site_coding == "lesion_dist_proximal") {
    adj_vars <- 
      adj_vars %>% 
      left_join(screening_data %>% 
                  inner_join(lesion_locs %>% 
                               select(most_serious_segment = segment, side_processed) %>% 
                               distinct() %>%
                               filter(side_processed == "Proximal") %>% 
                               na.omit(), by = join_by(most_serious_segment)) %>% 
                  select(most_serious_distancefromanus_processed, deltaker_id, final_result),
                by = join_by(deltaker_id)) %>% 
      select(-c(deltaker_id, final_result))
  }
  
  tmp_mod <- names(adj_vars %>% select(-sample_id)) %>% paste(collapse = "+")
  
  tmp_mod <- as.formula(paste("diversity ~", tmp_mod))
  
  
  alpha_diversity %>%
    pivot_longer(-c(sample_id, dataset), values_to = "diversity", names_to = "index") %>% 
    inner_join(adj_vars, by = "sample_id") %>% 
    group_by(index, dataset) %>% 
    group_split() %>% 
    lapply(function(ds) {
      
      lm(formula = tmp_mod, data = ds) %>% 
        summary() %>% 
        tidy(conf.int = TRUE) %>% 
        mutate(index = ds$index[1],
               dataset = ds$dataset[1])
    }) %>% 
    bind_rows() %>% 
    mutate(site_coding = site_coding)
}


## Presence/absence lesions by site
test_alpha_outcome_assoc_by_site_spec(site_coding = "lesion_loc_overall") %>% 
  ## Only assessed for those with distal lesions - not adjusted for presence of proximal lesions
  bind_rows(test_alpha_outcome_assoc_by_site_spec(site_coding = "lesion_loc_size_distal")) %>% 
  ## Only assessed for those with proximal lesions - not adjusted for presence of distal lesions
  bind_rows(test_alpha_outcome_assoc_by_site_spec(site_coding = "lesion_loc_size_proximal")) %>% 
  
  ## Only assessed for those with clinically relevant findings - and only for those with their most advanced finding distally
  bind_rows(test_alpha_outcome_assoc_by_site_spec(site_coding = "lesion_dist_distal")) %>%
  ## Only assessed for those with clinically relevant findings - and only for those with their most advanced finding proximally
  bind_rows(test_alpha_outcome_assoc_by_site_spec(site_coding = "lesion_dist_proximal")) %>%
  filter(str_detect(term, "(distal|proximal|most_serious|sum_lesion|sum_p|sum_d)")) #%>% 
  # write_tsv("results/tables/diversity_composition/alpha_site_specificity.tsv")





# Beta tests by lesion loc/burden -----------------------------------------

test_beta_outcome_assoc_by_site_spec <- function(dists, site_coding = "lesion_type_loc_dummy") {
  adj_vars <- 
    sample_data %>% 
    select(sample_id, deltaker_id) %>% 
    left_join(meta_dat %>% 
                select(deltaker_id, senter, age_invitation, 
                       kjonn, wcrf_index_main, Smoking, Utdanning,
                       colo_hemorrhoids, colo_IBD, antibiotics_reg_quest_comb), 
              by = "deltaker_id")
  
  if (site_coding == "lesion_type_loc_dummy") {
    adj_vars <- 
      adj_vars %>% 
      left_join(screening_data %>% 
                  select(contains("distal"), contains("proximal"), 
                         -contains("distance"), -contains("acn"), -contains("ssp_ssa"),
                         -starts_with("sum"),
                         deltaker_id), by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_loc_overall") {
    adj_vars <-
      adj_vars %>% 
      left_join(screening_data %>% 
                  select(contains("distal") & contains("acn"),
                         contains("proximal") & contains("acn"),
                         deltaker_id), by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_loc_size") {
    adj_vars <-
      adj_vars %>% 
      inner_join(screening_data %>% 
                   filter(!str_detect(final_result, "(1.|3.|4.)")) %>% 
                   select(contains("distal") & starts_with("sum"),
                          contains("proximal") & starts_with("sum"),
                          -contains("volume"), deltaker_id), 
                 by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_loc_size_distal") {
    adj_vars <-
      adj_vars %>% 
      inner_join(screening_data %>%
                   filter(sum_lesion_diameter_distal > 0) %>%
                   select(sum_lesion_diameter_distal,
                          deltaker_id),
                 by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_loc_size_proximal") {
    adj_vars <-
      adj_vars %>% 
      inner_join(screening_data %>%
                   filter(sum_lesion_diameter_proximal > 0) %>%
                   select(sum_lesion_diameter_proximal,
                          deltaker_id),
                 by = join_by(deltaker_id)) %>% 
      select(-deltaker_id)
  } else if (site_coding == "lesion_dist_distal") {
    adj_vars <- 
      adj_vars %>% 
      left_join(screening_data %>% 
                  inner_join(lesion_locs %>% 
                               select(most_serious_segment = segment, side_processed) %>% 
                               distinct() %>%
                               filter(side_processed == "Distal") %>% 
                               na.omit(), by = join_by(most_serious_segment)) %>% 
                  select(most_serious_distancefromanus_processed, deltaker_id, final_result),
                by = join_by(deltaker_id)) %>% 
      filter(!str_detect(final_result, "(1.|3.|4.)")) %>% 
      select(-c(deltaker_id, final_result))
  } else if (site_coding == "lesion_dist_proximal") {
    adj_vars <- 
      adj_vars %>% 
      left_join(screening_data %>% 
                  inner_join(lesion_locs %>% 
                               select(most_serious_segment = segment, side_processed) %>% 
                               distinct() %>%
                               filter(side_processed == "Proximal") %>% 
                               na.omit(), by = join_by(most_serious_segment)) %>% 
                  select(most_serious_distancefromanus_processed, deltaker_id, final_result),
                by = join_by(deltaker_id)) %>% 
      filter(!str_detect(final_result, "(1.|3.|4.)")) %>% 
      select(-c(deltaker_id, final_result))
  }
  tmp_sel <-
    dists %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column("sample_id") %>% 
    select(sample_id) %>% 
    left_join(adj_vars %>% 
                na.omit() %>% 
                mutate(sel = TRUE) %>% 
                select(sample_id, sel), by = join_by(sample_id)) %>% 
    mutate(sel = ifelse(is.na(sel), FALSE, sel)) %>% 
    pull(sel)
  
  tmp_meta <-   
    dists %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column("sample_id") %>% 
    select(sample_id) %>% 
    left_join(adj_vars, by = join_by(sample_id))
  
  print(paste0("Testing ", site_coding, ", n = ", sum(tmp_sel)))
  
  tmp <- vegan::adonis2(formula = as.dist(as.matrix(dists)[tmp_sel, tmp_sel]) ~ ., data = tmp_meta %>% filter(tmp_sel) %>% select(-sample_id), by = "margin")
  tmp2 <- adonis_OmegaSq(tmp) %>% as.data.frame() %>% rownames_to_column("term") %>% tibble() %>% bind_cols("R2" = tmp$R2)
  
  tmp2 %>% 
    mutate(n = sum(tmp_sel),
           site_coding = site_coding) %>% 
    filter(!term %in% c("Residual","Total"))
}

tmp_dists <-
  list("metaphlan" = mphlan_abundance,
       "mags" = MAGs_abundance,
       "keggs" = ko_hum) %>% 
  lapply(function(dataset) {
    dataset %>% 
      as.data.frame() %>% 
      column_to_rownames("sample_id") %>% 
      as.matrix() %>% 
      vegdist(method = "bray")
  })

beta_loc_tests <-
  test_beta_outcome_assoc_by_site_spec(tmp_dists[["metaphlan"]], site_coding = "lesion_loc_overall") %>% 
  bind_rows(test_beta_outcome_assoc_by_site_spec(tmp_dists[["metaphlan"]], site_coding = "lesion_loc_size_distal")) %>% 
  bind_rows(test_beta_outcome_assoc_by_site_spec(tmp_dists[["metaphlan"]], site_coding = "lesion_loc_size_proximal")) %>% 
  bind_rows(test_beta_outcome_assoc_by_site_spec(tmp_dists[["metaphlan"]], site_coding = "lesion_dist_distal")) %>% 
  bind_rows(test_beta_outcome_assoc_by_site_spec(tmp_dists[["metaphlan"]], site_coding = "lesion_dist_proximal"))

beta_loc_tests %>% 
  filter(str_detect(term, "(distal|proximal|most)")) #%>% 
  # write_tsv("results/tables/diversity_composition/beta_site_specificity.tsv")
