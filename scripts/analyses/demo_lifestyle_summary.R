

generate_demo_lifestyle_tables <- function() {
  
  add_na_repl <- function(x, repl) {
    ## Function to replace NA's in a factor with a replacement
    ifelse(is.na(x), repl, as.character(x)) %>% 
      factor(levels = unique(c(levels(x), repl)))
  }
  
  outc_table <- function(df, vars, outcome_var, gt = FALSE, header = "", missing_text = "Unknown", p = TRUE) {
    tab <- 
      df %>% 
      tbl_summary(include = c(all_of(vars)),
                  by = all_of(outcome_var), missing_text = missing_text,
                  percent = "row",
                  type = all_categorical() ~ "categorical") 
    if (p) tab <- tab %>% add_p()
    if (gt) tab <- tab %>% as_gt()
    if (header != "") tab <- tab %>% tab_header(header)
    return(tab)
  }
  
  
  rename_variables <- function(ds, renaming_df = variables) {
    ds %>% 
      rename_with(.fn = function(x) renaming_df %>% filter(var_id %in% x) %>% pull(var_name),
                  .cols = any_of(renaming_df$var_id))
  }
  
  select_ds <- function(incl_vars, ...) {
    
    sample_data %>% 
      select(sample_id, deltaker_id) %>% 
      left_join(meta_dat %>% 
                  rename_variables() %>% select(deltaker_id, `CRC-related findings`, `Colonoscopy result`, `Outcome`, 
                                                any_of(incl_vars)), by = "deltaker_id") %>% 
      rename_with(.fn = function(x) variables %>% filter(var_id %in% x) %>% pull(var_name),
                  .cols = any_of(variables$var_id)) %>% 
      mutate(across(.cols = any_of(incl_vars) & where(~ is.factor(.x) | is.character(.x)), 
                    .fns = function(x) rename_var_levels(x))) %>% 
      mutate(across(.cols = where(~ is.factor(.x) | is.character(.x)), .fns = function(x) fct_na_level_to_value(x, extra_levels = c("Missing", "Unknown"))))
    
  }
  
  
  # Overview of association with advanced findings --------------------------
  
  main_vars_v5 <- c("Age", "Sex", "Screening center", 
                    "FIT value",
                    "Programmatic screening round", "Previous negative FIT",
                    "Hemorrhoids", "Diverticulitis", "IBD",
                    "National background", "Education", "Employment status",
                    "WCRF", "Smoking", "Antibiotics use")
  
  gtsummary::theme_gtsummary_compact()
  main_tab <-
    main_vars_v5 %>% 
    select_ds() %>% 
    mutate(findings = ifelse(`CRC-related findings` == "positive", "CRC-related findings", "No related findings")) %>% 
    outc_table(vars = main_vars_v5, 
               outcome_var = "findings",
               missing_text = "No data/Missing/Unknown") %>% 
    as_gt() %>% 
    tab_header("Participant characteristics and CRC-related findings")
  gtsummary::reset_gtsummary_theme()
  
  gtsummary::theme_gtsummary_compact()
  main_tab_fin_res <-
    main_vars_v5 %>% 
    select_ds() %>% 
    outc_table(vars = main_vars_v5, 
               outcome_var = "Outcome",
               missing_text = "No data/Missing/Unknown",
               p = FALSE) %>% 
    as_gt() %>% 
    tab_header("Participant characteristics and clinicopathological diagnoses")
  gtsummary::reset_gtsummary_theme()
  
  main_tab %>% 
    gtsave(filename = "results/characteristics_by_crc_related_findings.html")
  
  main_tab_fin_res %>% 
    gtsave(filename = "results/characteristics_by_cp_diagnosis.html")
  
}

