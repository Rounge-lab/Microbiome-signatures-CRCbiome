dmms <- load_dmm(synthetic_data = TRUE)

i <- "final_result_cat_neg"
mod = "fully_adjusted_model"


## Conduct omnibus likelihood ratio tests for association with host variables:
## Two models: 1) simple host-association model, and 2) standard model

mnom_group_host_var_assoc <-
  dmms %>% 
  left_join(sample_data %>% select(sample_id, deltaker_id), by = "sample_id") %>% 
  left_join(meta_dat_cat %>% pivot_longer(-deltaker_id), by = "deltaker_id") %>% 
  select(-deltaker_id) %>% 
  # filter(name %in% "final_result") -> x
  group_by(name) %>% 
  group_split() %>% 
  lapply(function(x) {
    
    # lapply(c("fully_adjusted_model", "simplified_model", "sensitivity_no_crc"), function(mod) {
    lapply(c("host_assoc_model", "standard_model"), function(mod) {
      if (mod %in% "host_assoc_model") {
        adj_vars <- 
          sample_data %>% 
          select(sample_id, deltaker_id) %>% 
          left_join(meta_dat %>%
                      select(deltaker_id, senter, age_invitation, kjonn), by = "deltaker_id") %>%
          select(-c(deltaker_id, any_of(x$name[1])))
      }
      if (mod %in% "standard_model") {
        adj_vars <- 
          sample_data %>% 
          select(sample_id, deltaker_id) %>% 
          left_join(meta_dat %>% 
                      select(deltaker_id, senter, age_invitation, kjonn, 
                             wcrf_index_main, Smoking, Utdanning,
                             colo_hemorrhoids, colo_IBD, antibiotics_reg_quest_comb), by = "deltaker_id") %>%
          select(-c(deltaker_id, any_of(x$name[1])))
      }
      
      tmp_ref <- variables %>% filter(var_id %in% x$name[1]) %>% pull(ref)
      tmp_ref <- ifelse(is.na(tmp_ref), "low", tmp_ref)
      
      incl_mod <-
        x %>% 
        filter(!is.na(value)) %>% 
        inner_join(adj_vars, by = "sample_id") %>% 
        select(-c(name, sample_id)) %>% 
        mutate(value = fct_relevel(value, tmp_ref)) %>% 
        multinom(gr ~ ., data = .)  
      excl_mod <- 
        x %>% 
        filter(!is.na(value)) %>% 
        inner_join(adj_vars, by = "sample_id") %>% 
        select(-c(name, sample_id, value)) %>% 
        multinom(gr ~ ., data = .)  
      
      anova(excl_mod, incl_mod) %>% 
        tibble() %>% 
        mutate(var_id = x$name[1],
               model_complexity = mod)
    }) %>% 
      bind_rows()
  }) %>% 
  bind_rows()

mnom_group_host_var_assoc %>% 
  write_tsv("data/dmm/dmm_host_assoc_lrt.tsv")

## Plot
mnom_group_host_var_assoc %>% 
  filter(!is.na(`   Df`)) %>%
  mutate(model_complexity = factor(model_complexity, 
                                   levels = c("host_assoc_model", "standard_model"),
                                   labels = c("simplified model", "standard model"))) %>% 
  # mutate(model_complexity = case_when(str_detect(Model, "wcrf") ~ "simplified model",
  #                                     str_detect(Model, "Smoking") ~ "simplified model",
  #                                     TRUE ~ "base model")) %>% 
  arrange(desc(`Pr(Chi)`)) %>% 
  left_join(variables %>% select(var_id, var_name, lididem_crc), by = "var_id") %>% 
  filter(!var_id %in% c("final_result_cat_neg", "final_result", "detect_worthy_lesions")) %>% 
  View()
  filter(lididem_crc | 
           var_name %in% c("WCRF", "Screening center", "National background", "Employment status", "Education", "Antibiotics use")) %>% 
  mutate(var_name = factor(var_name, levels = unique(var_name[ model_complexity %in% "simplified model"]))) %>% 
  # filter(!str_detect(var_name, "gastrointestinal disorders")) %>% 
  ggplot(aes(x = -log10(`Pr(Chi)`), y = var_name, fill = model_complexity)) +
  geom_col(position = "dodge") +
  geom_vline(xintercept = -log10(0.05), color = "red", linetype = 2) +
  theme_bw() +
  scale_fill_manual(values = c(paletteer::paletteer_d(color_set)[c(2:3)], "gray")) +
  labs(y = "",
       fill = "") +
  facet_wrap(~model_complexity, nrow = 1)

dmms %>% 
  mutate(gr = fct_collapse(gr, test_level = "V3", ref_level = c("V1", "V2", "V4"))) %>% 
  count(gr)

## One v rest
dmm_mnom_outcome_assoc_one_v_rest <-
  lapply(paste0("V",1:4), function(dmm_test_lvl) {
    lapply(c("detect_worthy_lesions", "final_result_cat_neg"), function(i) {
      
      mod <- "standard_model"
      
      if (mod %in% "standard_model") {
        
        adj_vars <- 
          sample_data %>% 
          select(sample_id, deltaker_id) %>% 
          left_join(meta_dat %>% 
                      select(deltaker_id, senter, age_invitation, kjonn, 
                             wcrf_index_main, Smoking, Utdanning, 
                             colo_hemorrhoids, colo_IBD, antibiotics_reg_quest_comb,
                             outcome = any_of(i)), by = "deltaker_id") %>%
          select(-c(deltaker_id))
      }
      ref_lvls <- paste0("V", 1:4)
      ref_lvls <- ref_lvls[-which(ref_lvls %in% dmm_test_lvl)]
      
      tmp_mod <-
        dmms %>% 
        inner_join(adj_vars, by = "sample_id") %>%
        select(-sample_id) %>% 
        mutate(gr = fct_collapse(gr, test_level = dmm_test_lvl, ref_level = ref_lvls)) %>%
        mutate(gr = fct_relevel(gr, "ref_level")) %>%
        multinom(outcome ~ ., data = .) %>% 
        tidy(conf.int = TRUE, exponentiate = TRUE) %>% 
        filter(term != "(Intercept)") %>% 
        mutate(outcome = i,
               model = mod,
               dmm_test_lvl = dmm_test_lvl)
    }) %>% 
      bind_rows()
  }) %>% 
  bind_rows()

dmm_mnom_outcome_assoc_one_v_rest %>% 
  write_tsv("results/tables/dmm/dmm_mnom_outcome_assoc_one_v_rest.tsv")
