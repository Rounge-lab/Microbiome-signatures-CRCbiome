

diff_abund_meta_vars <- function() {
  # filter ko table ---------------------------------------------------------
  
  if (identical(rownames(ko_hum), rownames(mphlan_abundance))) {
    cor_matrix <- cor(ko_hum[,-1], mphlan_abundance[,-1])
    
    high_cor_paths <- colnames(ko_hum[,-1])[apply(cor_matrix, 1, function(x) any(x > 0.5))]
    
    ko_hum_filtered <- 
      ko_hum %>% 
      select(-all_of(high_cor_paths))
  } else {
    stop("make sure data is aligned")
  }
  
  # run maaslin -------------------------------------------------------------
  
  sink("data/diff_abund/Maaslin2/meta_vars/log.log")
  
  diff_abund_results <-
    lapply(c("mphlan", "MAGs", "ko_hum", "genera"), function(dataset) {
      if (dataset == "mphlan") tmp_abund <- mphlan_abundance
      if (dataset == "MAGs") tmp_abund <- MAGs_abundance
      if (dataset == "ko_hum") tmp_abund <- ko_hum_filtered
      
      tmp_abund <- tmp_abund %>% mutate(sample_id = str_replace(sample_id, "-", "_"))
      
      if (!dataset %in% list.dirs("data/diff_abund/Maaslin2/meta_vars/", full.names = FALSE)) dir.create(paste0("data/diff_abund/Maaslin2/meta_vars/", dataset))
      
      meta_dat %>% 
        inner_join(sample_data %>% select(sample_id, deltaker_id)) %>% 
        select(-deltaker_id) %>% 
        mutate(across(where(is.double), .fns = function(x) cat_func(x) %>% str_replace("negative", "a_low") %>% str_replace("positive", "b_high"))) %>% 
        pivot_longer(-sample_id) %>% 
        group_by(name) %>% 
        group_split() %>% 
        lapply(function(x) {
          
          tmp_meta <- x
          
          tmp_ref <- variables %>% mutate(ref = ifelse(is.na(ref), "a_low", ref)) %>% filter(var_id %in% x$name[1]) %>% pull(ref)
          
          adj_vars <- c("value", "kjonn","age_invitation","senter")
          
          x <- 
            tmp_meta %>% 
            left_join(meta_dat %>% 
                        inner_join(sample_data %>% select(sample_id, deltaker_id), by = "deltaker_id") %>% 
                        select(sample_id, any_of(adj_vars)) %>% 
                        select(-any_of(x$name[1])), by = "sample_id")
          
          mod <- paste0(adj_vars[!adj_vars %in% x$name[1]], collapse = ",")
          refs <- tibble(variable = c("value", "kjonn", "senter"), 
                         ref = c(tmp_ref, levels(meta_dat$kjonn)[1], levels(meta_dat$senter)[1])) %>% 
            filter(!variable %in% x$name[1]) %>%
            mutate(f = paste(variable, ref, sep = ",")) %>% 
            pull(f) %>% 
            paste(collapse = ";")
          
          
          tmp_maaslin <-
            Maaslin2::Maaslin2(
              input_data = x %>% 
                select(sample_id) %>% 
                left_join(tmp_abund) %>% 
                as.data.frame() %>% 
                column_to_rownames("sample_id"),
              input_metadata = x %>% 
                column_to_rownames("sample_id"),
              min_prevalence = 0.1,
              normalization = "TSS",
              transform = "LOG",
              analysis_method = "LM",
              output = paste0("data/diff_abund/Maaslin2/meta_vars/", dataset, "/", x$name[1]),
              ## Add adjustment: Age + sex
              fixed_effects = mod,
              reference = refs,
              plot_heatmap = FALSE,
              plot_scatter = FALSE
            )
          
          tmp_maaslin$results %>% 
            tibble() %>% 
            mutate(variable = x$name[1],
                   dataset = dataset) %>% 
            filter(metadata %in% "value")
          
        }) %>% 
        bind_rows()
    }) %>% 
    bind_rows()
  
  sink()
    
  diff_abund_results %>%
    write_tsv("data/diff_abund/Maaslin2/meta_vars/maaslin2_meta_vars.tsv")
  
}

diff_abund_meta_vars()

rm(diff_abund_meta_vars)
