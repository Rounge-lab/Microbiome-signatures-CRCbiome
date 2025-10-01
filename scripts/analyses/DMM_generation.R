


run_dmm <- function(dataset = "crcbiome", k = 1:10) {
  
  if (dataset == "crcbiome") {
    crcbiome_filtering <- 
      mphlan_abundance %>% 
      # slice_sample(n = 150) %>% 
      as.data.frame() %>% 
      column_to_rownames("sample_id") %>% 
      select(which(apply(. > 0.1, 2, mean) > 0.1)) %>% 
      tibble() %>% 
      names()
    
    ## CRCbiome - will set the same parameters
    
    crcbiome_dmms <- 
      lapply(k, function(x) {
        # lapply(4, function(x) {
        print(x)
        set.seed(1)
        mphlan_abundance %>% 
          as.data.frame() %>% 
          column_to_rownames("sample_id") %>% 
          select(which(apply(. > 0.1, 2, mean) > 0.1)) %>% 
          as.matrix() %>% 
          DirichletMultinomial::dmn(k = x, verbose = TRUE)
      })
    
    crcbiome_dmms %>% 
      write_rds("data/dmm/crcbiome_dmm.rds")
   
  } else if (dataset == "lee") {
    lee_meta <- read_tsv("data/input_processed/lee_metadata.tsv", col_types = cols())
    lee_sgbs <- read_rds("data/input_processed/lee_sgbs.Rds")
    lee_sgb_names <- read_tsv("data/input_processed/lee_sgb_names.tsv", col_types = cols())
    
    
    lee_filtering <- 
      lee_sgbs %>% 
      # slice_sample(n = 150) %>% 
      as.data.frame() %>% 
      column_to_rownames("sample_id") %>% 
      select(which(apply(. > 0.1, 2, mean) > 0.1)) %>% 
      tibble() %>% 
      names()
    
    lee_dmms <- 
      lapply(k, function(x) {
        print(x)
        set.seed(1)
        lee_sgbs %>% 
          # slice_sample(n = 150) %>% 
          as.data.frame() %>% 
          column_to_rownames("sample_id") %>% 
          select(which(apply(. > 0.1, 2, mean) > 0.1)) %>% 
          as.matrix() %>% 
          DirichletMultinomial::dmn(k = x, verbose = TRUE)
      })
    
    lee_dmms %>% 
      write_rds("data/dmm/lee_dmm.rds")
  }
  
}

evaluate_dmms <- function(dmm_file, dataset = "crcbiome") {
  tmp <- read_rds(dmm_file) %>% 
    sapply(DirichletMultinomial::laplace)
  pdf(paste0("data/dmm/", dataset, "_laplace.pdf"))
  plot(tmp, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit", main = dataset)
  dev.off()
}

define_dmm_group <- function(dmm_file, k, dataset = "crcbiome") {
  read_rds(dmm_file)[[k]]@group %>% 
    as.data.frame() %>% 
    mutate(gr = paste0("V",seq(ncol(.)))[ apply(., 1, which.max)]) %>% 
    rownames_to_column("sample_id") %>% 
    tibble() %>% 
    select(-starts_with("V", ignore.case = FALSE)) 
}

get_posterior_mean_diff <- function(dmm_file = "data/dmm/crcbiome_dmm.rds", k = 4, dataset = "crcbiome") {
  
  dmm <- read_rds(dmm_file)[[k]]
  base_dmm <- read_rds(dmm_file)[[1]]
  
  posterior_mean_diff <-
    base_dmm %>% 
    DirichletMultinomial::fitted(scale = TRUE) %>%
    as.data.frame() %>% 
    rownames_to_column("sgb") %>% 
    tibble() %>% 
    rename(single_comp = 2) %>% 
    left_join(dmm %>% 
                DirichletMultinomial::fitted(scale = TRUE) %>%
                as.data.frame() %>% 
                rownames_to_column("sgb") %>% 
                tibble(), by = "sgb") %>% 
    mutate(diff = apply(abs(single_comp - .[,paste0("V",seq(k))]), 1, sum)) %>% 
    arrange(desc(diff)) %>% 
    mutate(cdiff = cumsum(diff/sum(diff))) %>% 
    mutate(max_g = paste0("V",seq(k))[ apply(abs(single_comp-.[,paste0("V",seq(k))]), 1, which.max)]) %>% 
    pivot_longer(all_of(paste0("V", seq(k))), names_to = "component", values_to = "value") %>%
    mutate(rel_diff = log10(value/single_comp)) %>%
    pivot_wider(names_from = component, values_from = c(value, rel_diff)) %>%
    left_join(mp4_taxonomy %>% select(sgb, species)) %>%
    mutate(alt_max_g = paste0("V",seq(k))[ apply(abs(single_comp-.[,paste0("rel_diff_V",seq(k))]), 1, which.max)]) %>%
    rename_with(.fn = function(x) str_remove(x, "value_"), .cols = starts_with("value"))
  
  posterior_mean_diff %>% 
    write_tsv("results/tables/dmm/posterior_mean_diff.tsv")
}

run_dmm()
# run_dmm(dataset = "lee")
evaluate_dmms("data/dmm/crcbiome_dmm.rds", dataset = "crcbiome")
# evaluate_dmms("data/dmm/lee_dmm.rds", dataset = "lee")
## Define groups for minimum model error
define_dmm_group("data/dmm/crcbiome_dmm.rds", k = 4, dataset = "crcbiome") %>%
  write_tsv("data/dmm/crcbiome_dmm.tsv")
get_posterior_mean_diff()


