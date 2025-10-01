

load_data <- function(synthetic_data = FALSE) {
  if (synthetic_data) {
    source("scripts/analyses/load_synthetic_data.R")
  } else {
    source("scripts/analyses/load_data.R")  
  }
}

load_dmm <- function(synthetic_data = FALSE) {
  if (synthetic_data) {
    read_tsv("data/synthetic/crcbiome_dmm.tsv", col_types = cols())
  } else {
    read_tsv("data/dmm/crcbiome_dmm.tsv", col_types = cols())
  }
}

load_colibactin_data <- function(synthetic_data = FALSE) {
  if (synthetic_data) {
    read_tsv("data/synthetic/pks_detection.tsv", col_types = cols())
  } else {
    read_tsv("data/input_processed/pks_detection.tsv", col_types = cols())
  }
}


## Define set of var levels to consider
define_prev_var_levels <- function(dat = meta_dat_cat, vars = variables, min_n = 10, var_name_and_id = "both", name_update = TRUE) {
  
  tmp_name <- 
    dat %>% 
    pivot_longer(-deltaker_id) %>% 
    count(name, value) %>% 
    (function(x) {
      if (name_update) x %>% mutate(value = rename_var_levels(value))
      else x
    }) %>% 
    filter(!is.na(value)) %>% 
    left_join(vars %>% select(name = var_id, var_name), by = "name") %>% 
    mutate(var = paste(var_name, value)) %>% 
    filter(n >= min_n) %>% 
    pull(var) 
  tmp_id <-
    dat %>% 
    pivot_longer(-deltaker_id) %>% 
    count(name, value) %>% 
    (function(x) {
      if (name_update) x %>% mutate(value = rename_var_levels(value))
      else x
    }) %>% 
    filter(!is.na(value)) %>% 
    mutate(var = paste(name, value)) %>% 
    filter(n >= min_n) %>% 
    pull(var)
  
  if (var_name_and_id %in% "both") {
    c(tmp_name, tmp_id)  
  } else if (var_name_and_id %in% "name") {
    tmp_name
  } else if (var_name_and_id %in% "id") {
    tmp_id
  }
} 

## Function to make a tibble that can be used to create a PCoA plot
dist_to_PCoA <- function(distance_matrix, group_var = "group", PCoA_dims = c(1,2), return_centroids = TRUE) {
  if (return_centroids) {
    
    pcoa_obj <- betadisper(d = distance_matrix, group = group_var)
    
    tmp <-
      pcoa_obj$vectors %>% 
      data.frame() %>% 
      select(all_of(PCoA_dims)) %>% 
      rownames_to_column("sample_id") %>%
      tibble() %>% 
      bind_cols(group = pcoa_obj$group) %>% 
      pivot_longer(-c(sample_id, group), names_to = "PCoA", values_to = "PCo_value") %>% 
      left_join(pcoa_obj$centroids %>% 
                data.frame() %>% 
                select(PCoA_dims) %>% 
                rownames_to_column("group") %>% 
                tibble() %>% 
                pivot_longer(-c(group), names_to = "PCoA", values_to = "centroid")) %>% 
      left_join(pcoa_obj$eig %>% 
                  enframe(name = "PCoA") %>% 
                  mutate(rel_eig = value/sum(value)) %>% 
                  select(PCoA, rel_eig)) %>% 
      pivot_longer(c(PCo_value, centroid), names_to = "var_type")
  } else {
    pcoa_obj <- 
      distance_matrix %>% 
      pcoa()
    tmp <-
      pcoa_obj$vectors %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("sample_id") %>% 
      tibble() %>% 
      select(sample_id, PCoA_dims+1) %>% 
      pivot_longer(-sample_id, values_to = "PCo_value", names_to = "PCoA")
  }
  tmp
}

## Function to plot two PCoA dimensions
plot_pcoa <- function(x, dim_1 = "PCoA1", dim_2 = "PCoA2", type = "star", with_bars = FALSE) {
  tmp_labs <- x %>% 
    filter(PCoA %in% c(dim_1, dim_2)) %>% 
    mutate(lab = paste(str_replace(PCoA, "PCoA", "PCo"), " (", round(rel_eig*100, digits = 1), "%)", sep = "")) %>% 
    group_by(lab) %>% 
    slice_head(n = 1) %>% 
    ungroup() %>% 
    pull(lab)
  
  x <- x %>% 
    select(-rel_eig) %>% 
    pivot_wider(names_from = PCoA, values_from = value) %>%
    rename_with(.fn = function(x) ifelse(x == dim_1, "dim1", x)) %>% 
    rename_with(.fn = function(x) ifelse(x == dim_2, "dim2", x))
  
  if (type %in% "star") {
    tmp_plot <-
      x %>% 
      ggplot(aes(x = dim1, y = dim2, color = group, group = sample_id)) +
      geom_line(alpha = .1, show.legend = FALSE) +
      geom_point(data = x %>% filter(var_type %in% "centroid"), aes(x = dim1, y = dim2, fill = group), color = "black", shape = 21, size = 3) +
      stat_ellipse(data = x %>% filter(var_type %in% "PCo_value"), aes(x = dim1, y = dim2, color = group, group = group), level = .5, show.legend = FALSE)
      
  } else {
    tmp_plot <-
      x %>% 
      ggplot(aes(x = dim1, y = dim2, color = group, group = sample_id)) +
      geom_point(alpha = .2, show.legend = FALSE) +
      geom_point(data = x %>% filter(var_type %in% "centroid"), aes(x = dim1, y = dim2, fill = group), color = "black", shape = 21, size = 3) +
      stat_ellipse(data = x %>% filter(var_type %in% "PCo_value"), aes(x = dim1, y = dim2, color = group, group = group), level = .5, show.legend = FALSE)
  }
  
  tmp_plot <- 
    tmp_plot +
    theme_bw() +
    labs(color = "", fill = "", x = dim_1, y = dim_2)  +
    labs(x = tmp_labs[1],
         y = tmp_labs[2]) 
  
  if (with_bars) {
    tmp_plot_1 <-
      x %>% 
      filter(var_type %in% "PCo_value") %>% 
      ggplot(aes(x = dim1, y = 1, fill = group)) +
      geom_boxplot(outlier.shape = 20, outlier.size = 0.2) +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none") +
      labs(y = "", x = "")
    
    tmp_plot_2 <-
      x %>% 
      filter(var_type %in% "PCo_value") %>% 
      ggplot(aes(x = 1, y = dim2, fill = group)) +
      geom_boxplot(outlier.shape = 20, outlier.size = 0.2) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none") +
      labs(y = "", x = "")
    
    list(tmp_plot_1, NULL, tmp_plot, tmp_plot_2)
    
  } else {
    tmp_plot +
      coord_fixed()
  }
}

## micEco function to obtain omega^2 - effect size - for permanova

# source("workflow/scripts/micEco/adonis_Omega.R")


## Function to categorize continuous variables into tertiles - except in cases where > 33% of cases are 0, when tertiles are based on the positive cases.
cat_func <- function(x, labs = c("negative", "mid", "positive")) {
  if(is.numeric(x)) {
    # print("cat")
    ## median for those with more than 0 - if these comprise more than 33%
    if (sum(x %in% 0) >= 0.33*sum(!is.na(x))) {
      cut(x, 
          breaks = c(-Inf, 0, quantile(x[ x > 0], probs = seq(0, 1, length.out = 3), na.rm = TRUE, names = FALSE)[2], Inf), 
          labels = labs)  
    } else {
      ## Otherwise, tertiles for whole range
      cut(x, 
          breaks = c(-Inf, quantile(x, probs = seq(0, 1, length.out = 4), na.rm = TRUE, names = FALSE)[2:3], Inf), 
          labels = labs)  
    }
  } else {
    x
  }
}

## Order column names of matrix/dataframe by hierarchical clustering
order_by_cluster_columns <- function(x) {
  
  tmp_0 <- names(x)[which(apply(x == 0, 2, all))]
  x <- x[,!apply(x == 0, 2, all)]
  
  tmp <- x %>% scale() %>% t() %>% dist() %>% hclust()
  c(names(x)[tmp$order], tmp_0) 
  
}

## Order column names of matrix/dataframe by hierarchical clustering
order_by_cluster_variable <- function(x, var_tobe_clustered, var_to_cluster_on, value_to_cluster_on) {
  
  x %>% 
    select(all_of(c(var_tobe_clustered, var_to_cluster_on, value_to_cluster_on))) %>% 
    rename(x_var = all_of(var_tobe_clustered),
           y_var = all_of(var_to_cluster_on),
           n_var = all_of(value_to_cluster_on)) %>% 
    pivot_wider(names_from = x_var, values_from = n_var, values_fill = 0) %>% 
    as.data.frame() %>% 
    column_to_rownames("y_var") %>% 
    order_by_cluster_columns()
}


## Create sensitivity, specificity plot by threshold with confidence intervals
roc_plot_obj <- function(dat, 
                         response = "detect_worthy_lesions", 
                         predictor = ".pred_positive",
                         resp_levels = c("negative", "positive")) {
  
  dat <- dat %>% 
    rename(resp = all_of(response), 
           pred = all_of(predictor))
  
  tmp_roc <- dat %>% 
    roc(response = resp, 
        predictor = pred, 
        levels = resp_levels, 
        direction = "<")
  
  tmp_ci <- ci.se(tmp_roc, specificities = tmp_roc$specificities, boot.n = 1000)
  
  
  tibble(sensitivities = tmp_roc$sensitivities,
         specificities = tmp_roc$specificities,
         thresholds = tmp_roc$thresholds,
         var_cat = "roc") %>%
    bind_rows(tibble(specificities = as.numeric(row.names(tmp_ci)),
                     lower_ci = tmp_ci[,1],
                     upper_ci = tmp_ci[,3]) %>%
                unique() %>% 
                pivot_longer(-specificities, 
                             names_to = "var_cat",
                             values_to = "sensitivities"))
  
}


## Get some theme color
col_theme <- function(theme_col, n, preview = FALSE) {
  tmp <- paletteer::paletteer_d(theme_col)[1:n]  
  if (!preview) {
    tmp %>% 
      as.character() %>% 
      return()
  } else {
    tmp %>% 
      return()
  }
}

if (! "color_set" %in% ls()) color_set <- "ggsci::nrc_npg"

color_assignments <- 
  tibble(variable = 
           c("1. Negative", 
             "3. Non neoplastic findings", 
             "4. Other lesions", 
             "5b. >= 3 Non-advanced adenomas", 
             "5c. Advanced adenoma", 
             "5d. Advanced serrated", 
             "6. Cancer"),
         color = c(paletteer::paletteer_d(color_set)[c(9,5,10,4:1)])) %>% 
  bind_rows(tibble(variable = c("Not clinically relevant", 
                                "Not CRC-related",
                                "Clinically relevant",
                                "Clinically relevant findings",
                                "CRC-related findings",
                                "CRC-related"),
                   color = c(rep("gray", 2), paletteer::paletteer_d(color_set)[c(6,6,6,6)]))) %>% 
  bind_rows(tibble(variable = c("Negative", 
                                ">= 3 Non-advanced adenomas", 
                                "Advanced adenoma", 
                                "Advanced serrated", 
                                "Cancer"),
                   color = c("gray", paletteer::paletteer_d(color_set)[c(5,3:1)]))) %>% 
  bind_rows(tibble(variable = c("microbial model alone", "microbial model + FIT", 
                                "FIT alone", "microbial model + FIT, sex, age, and WCRF score", 
                                "FIT, sex, age, and WCRF score"),
                   color = col_theme("fishualize::Etheostoma_spectabile", 5))) %>% 
  bind_rows(tibble(variable = c("P. copri enriched", "M. smithii enriched",
                                "Clostridia enriched", "Bacteroides enriched"),
                   color = col_theme("suffrager::CarolMan", 4))) %>% 
  deframe()
# paletteer::palettes_d_names %>% view()

## Function for renaming of variable levels for presentation
rename_var_levels <- function(var) {
  var_names <- c("Moss" = "Center 1",
                 "Bærum" = "Center 2",
                 "Native" = "Norwegian",
                 "Non-native" = "Non-Norwegian",
                 "detect_worthy_lesions" = "Clinically relevant findings",
                 "final_result" = "Colonoscopy result",
                 "final_result_cat_neg" = "Outcome category",
                 "V1" = "Clostridia enriched",
                 "V2" = "M. smithii enriched",
                 "V3" = "Bacteroides enriched",
                 "V4" = "P. copri enriched")
  var %>% 
    as.character() %>% 
    enframe(value = "orig") %>% 
    left_join(var_names %>% 
                enframe(name = "orig",
                        value = "new"), by = "orig") %>% 
    mutate(new = ifelse(is.na(new), orig, new)) %>% 
    pull(new)
}

## Function to annotate and categorize models
make_model_summary <- function(prob_list) {
  prob_list %>%
    pull(model) %>%
    unique() %>%
    enframe(value = "model") %>%
    select(-name) %>%
    mutate(data_type = case_when(grepl("mag", model) ~ "mag-based",
                                 grepl("hum", model) | grepl("mphl", model) ~ "read-based",
                                 TRUE ~ "other"),
           dataset_type = case_when(grepl("mag", model) ~ "mags",
                                 grepl("hum", model) ~ "genes",
                                 grepl("mphl", model) ~ "sgbs",
                                 TRUE ~ "other"),
           stratification = factor(case_when(grepl("FIT-high", model) ~ "FIT-high",
                                             grepl("FIT-low", model) ~ "FIT-low",
                                             grepl("WCRF-high", model) ~ "WCRF-high",
                                             grepl("WCRF-low", model) ~ "WCRF-low",
                                             grepl("-men", model) ~ "men",
                                             grepl("women", model) ~ "women",
                                             grepl("old", model) ~ "old",
                                             grepl("young", model) ~ "young",
                                             grepl("neg-v-crc$", model) ~ "neg-v-crc",
                                             grepl("neg-v-aa$", model) ~ "neg-v-aa",
                                             grepl("neg-v-aa-crc$", model) ~ "neg-v-aa-crc",
                                             grepl("neg-v-aa-serr-crc$", model) ~ "neg-v-aa-as-crc",
                                             grepl("neg-v-naa-aa-serr-crc$", model) ~ "neg-v-naa-aa-as-crc",
                                             grepl("neg-nnl-v-aa-serr-crc$", model) ~ "neg-nnl-v-aa-as-crc",
                                             grepl("neg-v-serr$", model) ~ "neg-v-as",
                                             grepl("proximal", model) ~ "proximal",
                                             grepl("distal$", model) ~ "distal",
                                             grepl("_mags$", model) ~ "none",
                                             grepl("_mphl$", model) ~ "none"),
                                   levels = c("old", "young", "men", "women", "WCRF-high", "WCRF-low", "proximal", "distal", "FIT-high", "FIT-low",
                                              "neg-v-aa-as-crc", "neg-v-naa-aa-as-crc", "neg-nnl-v-aa-as-crc", "neg-v-as", "neg-v-aa", "neg-v-aa-crc", "neg-v-crc", "none")),
           strat_cat = case_when(grepl("FIT-", stratification) ~ "FIT",
                                 grepl("WCRF-", stratification) ~ "WCRF",
                                 stratification %in% c("men", "women") ~ "sex",
                                 stratification %in% c("old", "young") ~ "age",
                                 stratification %in% c("distal", "proximal") ~ "localization",
                                 grepl("-v-", stratification) ~ "outcome",
                                 grepl("_mags$", model) ~ "none",
                                 grepl("_mphl$", model) ~ "none") %>% 
             factor(levels = c("FIT", "WCRF", "sex", "age", "localization", "outcome", "none")),
           dataset = factor(case_when(grepl("_mags$", model) | grepl("_mphl$", model) ~ "taxa",
                                      grepl("_mag-kegg$", model) | grepl("_hum-kegg$", model) ~ "genes",
                                      grepl("_lididem$", model) ~ "life, demo",
                                      grepl("_mags-mag-kegg$", model) | grepl("_mphl-hum-kegg$", model) ~ "taxa, genes",
                                      grepl("_mags-lididem$", model) | grepl("_mphl-lididem$", model) ~ "taxa, life, demo",
                                      grepl("_mag-kegg-lididem$", model) | grepl("_hum-kegg-lididem$", model) ~ "genes, life, demo",
                                      grepl("_mags-mag-kegg-lididem$", model) | grepl("_mphl-hum-kegg-lididem$", model) ~ "taxa, genes, life, demo",
                                      grepl("_mags-fit", model) | grepl("_mphl-fit", model) ~ "taxa, fit",
                                      grepl("_mags-lididem-fit", model) | grepl("_mphl-lididem-fit", model) ~ "taxa, life, demo, fit",
                                      grepl("_mags-wcrf-demo-fit", model) | grepl("_mphl-wcrf-demo-fit", model) ~ "taxa, wcrf, demo, fit"),
                            levels = rev(c("taxa", "genes", "life, demo",
                                           "taxa, genes", "taxa, life, demo", 
                                           "genes, life, demo", "taxa, genes, life, demo",
                                           "taxa, fit", "taxa, life, demo, fit", "taxa, wcrf, demo, fit"))),
           model_name = case_when(is.na(dataset) ~ as.character(stratification),
                                  TRUE ~ as.character(dataset)),
           algorithm = str_extract(model, "^.*(?=_)") %>% 
             factor(levels = c("lasso", "nnet", "rf", "svm", "xgb"), 
                    labels = c("lasso", "neural net", "random forest", "support vector machine", "extreme gradient boosting")))
}


summarize_model_stats <- function(grouped_data, stats = "auroc", output_wide = TRUE, join_by_var = "mod") {
  lapply(stats, function(stat) summarize_metric(grouped_data, stat, names_wide = output_wide)) %>% 
    (function(x) {
      if (output_wide) {
        x %>% reduce(left_join, by = join_by_var)
      } else {
        x %>% bind_rows()
      }
    })
}

## Summarize model metric across n times repeated k-folds
summarize_metric <- function(grouped_data, var = "auroc", aggregation = "n_by_k", names_wide = TRUE) {
  
  grouped_data <-
    grouped_data %>% 
    rename(metric_var = !!var)
  
  if (aggregation == "n_by_k") {
    grouped_data %>% 
      summarize(mean = mean(metric_var, na.rm = TRUE),
                sd = sd(metric_var, na.rm = TRUE),
                ci_low = mean - (sd(metric_var, na.rm = TRUE)*1.96/sqrt(n())),
                ci_high = mean + (sd(metric_var, na.rm = TRUE)*1.96/sqrt(n())),
                .groups = "drop") %>%
      (function(x) {
        if (names_wide) {
          x %>% rename_with(.fn = function(x) paste0(var, "_", x), .cols = matches("(mean|sd|ci_low|ci_high)"))
        } else {
          x %>% mutate(measure = var)
        }
      })
      
  } else if (aggregation == "by_rep") {
    grouped_data %>% 
      nest() %>% 
      mutate(summary = map(.x = data, .f = function(x) {
        x %>% 
          group_by(rep) %>% 
          summarize(mean = mean(metric_var, na.rm = TRUE),
                    .groups = "drop") %>% 
          summarize(m_mean = mean(mean),
                    sem = sd(mean),
                    ci_low = m_mean - sd(mean)*1.96/sqrt(n()),
                    ci_high = m_mean + sd(mean)*1.96/sqrt(n())) %>% 
          (function(x) {
            if (names_wide) {
              x %>% rename_with(.fn = function(x) paste0(var, "_", x), .cols = matches("(m_mean|sem|ci_low|ci_high)"))
            } else {
              x %>% mutate(measure = var)
            }
          })
      })) %>% 
      select(-data) %>% 
      unnest(summary)
  }
  
}
