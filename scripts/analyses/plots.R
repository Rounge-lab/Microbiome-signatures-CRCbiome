
plot_microbiome_overview <- function() {
  
  
  ds <- "metaphlan"
  
  tmp_var_levels <-
    define_prev_var_levels(var_name_and_id = "both")
  
  tmp_abund <- mphlan_abundance
  
  lab_tab <-
    mp4_taxonomy %>%
    select(feature = sgb, lab = species)
  
  lab_tab <-
    tmp_abund %>% 
    pivot_longer(-sample_id, values_to = "abundance", names_to = "feature") %>% 
    left_join(mp4_taxonomy %>% select(phylum, feature = sgb), by = "feature") %>% 
    group_by(sample_id, phylum) %>% 
    summarize(abundance = sum(abundance), .groups = "drop") %>% 
    filter(abundance > 0) %>% 
    mutate(lab = fct_infreq(phylum) %>% fct_lump_n(10) %>% fct_rev()) %>% 
    distinct(phylum, lab) %>% 
    left_join(mp4_taxonomy %>% select(feature = sgb, phylum), by = "phylum")
  
  prev_abundance_by_outcome <- 
    tmp_abund %>% 
    pivot_longer(-sample_id, values_to = "abundance", names_to = "feature") %>% 
    left_join(lab_tab, by = "feature") %>%
    group_by(lab, sample_id) %>% 
    summarize(abundance = sum(abundance), .groups = "drop") %>% 
    group_by(lab) %>% 
    summarize(prevalence_all = sum(abundance > 0)/n(),
              mean_abundance_all = mean(abundance[abundance > 0]), 
              .groups = "drop") %>% 
    ungroup() %>% 
    left_join(tmp_abund %>% 
                pivot_longer(-sample_id, values_to = "abund", names_to = "feature") %>% 
                left_join(lab_tab, by = "feature") %>% 
                group_by(lab, sample_id) %>% 
                summarize(abundance = sum(abund),
                          n_spec_pres = sum(abund>0), .groups = "drop") %>% 
                left_join(sample_data %>% 
                            left_join(screening_data %>% 
                                        select(deltaker_id, detect_worthy_lesions, final_result), 
                                      by = "deltaker_id") %>% 
                            select(sample_id, detect_worthy_lesions, final_result), 
                          by = "sample_id"), by = "lab") 
  
  prev_summary_by_outcome <-
    prev_abundance_by_outcome %>% 
    group_by(lab, detect_worthy_lesions, prevalence_all, mean_abundance_all) %>% 
    summarize(mean_abunance = mean(abundance[ abundance > 0]),
              prevalence = sum(abundance > 0)/n(), .groups = "drop")
  
  prev_plot <-
    prev_summary_by_outcome %>% 
    mutate(Findings = case_when(detect_worthy_lesions == "negative" ~ "Not clinically relevant",
                                detect_worthy_lesions == "positive" ~ "Clinically relevant")) %>% 
    ggplot(aes(x = prevalence, y = lab, fill = Findings)) +
    geom_col(position = position_dodge()) +
    labs(x = "Prevalence", y = "") +
    scale_fill_manual(values = color_assignments) +
    theme_bw()
  
  prev_abundance_plot <-
    prev_abundance_by_outcome %>% 
    group_by(lab) %>% 
    mutate(abundance = abundance / 100) %>% 
    mutate(ab_adj = abundance + min(abundance[ abundance > 0])/2) %>% 
    ungroup() %>% 
    filter(abundance > 0) %>% 
    mutate(Findings = case_when(detect_worthy_lesions == "negative" ~ "Not clinically relevant",
                                detect_worthy_lesions == "positive" ~ "Clinically relevant")) %>% 
    ggplot(aes(x = ab_adj, y = lab, fill = Findings)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
    scale_x_continuous(trans = "log10") +
    theme_bw() +
    labs(y = "", x = "Relative abundance*") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_fill_manual(values = color_assignments)
  
  species_per_phylum_plot <-
    prev_abundance_by_outcome %>% 
    mutate(Findings = case_when(detect_worthy_lesions == "negative" ~ "Not clinically relevant",
                                detect_worthy_lesions == "positive" ~ "Clinically relevant")) %>% 
    ggplot(aes(x = n_spec_pres, y = lab, fill = Findings)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
    scale_x_continuous(trans = "log10") +
    theme_bw() +
    labs(y = "", x = "Number of observed species*") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_fill_manual(values = color_assignments)
  
  prev_ab_plot <- ggpubr::ggarrange(prev_plot + theme(text = element_text(size = 6)), 
                                    prev_abundance_plot + theme(text = element_text(size = 6)),
                                    species_per_phylum_plot + theme(text = element_text(size = 6)),
                                    widths = c(1.4, 1, 1),
                                    common.legend = TRUE, legend = "bottom", labels = c("a", "b", "c"), 
                                    align = "h", nrow = 1)
  
  
  alpha_tests_cat <- 
    read_tsv("results/tables/alpha_host_vars_cat.tsv", col_types = cols()) %>% 
    filter(dataset %in% ds) %>% 
    filter(variable == test_var) %>% 
    left_join(variables %>% select(variable = var_id, var_name, variable_category = dataset, lididem_crc), by = "variable") %>% 
    mutate(lvl = rename_var_levels(lvl)) %>% 
    mutate(var_level = paste(var_name, lvl)) %>% 
    mutate(index = factor(index, levels = c("observed", "shannon", "invsimpson"))) 
  
  tmp_alpha_tests_cat <-
    alpha_tests_cat %>% 
    filter(var_level %in% tmp_var_levels) %>% 
    group_by(test_var) %>% 
    mutate(m_est = mean(estimate[ index %in% "observed" & !lvl %in% c("Unknown", "Missing")])) %>% 
    ungroup() %>% 
    arrange(m_est, index, estimate) %>% 
    mutate(var_level = factor(var_level, levels = var_level[ index %in% "observed"])) #%>% 
  
  alpha_assoc_plot <-
    tmp_alpha_tests_cat %>% 
    filter(!test_var %in% c("final_result", "detect_worthy_lesions")) %>% 
    filter(lididem_crc | 
             var_name %in% c("WCRF", "Screening center", "National background", 
                             "Employment status", "Education", "Antibiotics use",
                             "Diverticulitis", "Hemorrhoids", "IBD"),
           lvl != "Missing") %>%
    ggplot(aes(x = estimate, y = var_level)) +
    geom_linerange(aes(xmin = conf.low, xmax = conf.high), color = scales::muted("red")) +
    geom_point(color = "black", aes(shape = p.value < 0.05)) +
    scale_shape_manual(values = c(1,20)) +
    geom_vline(xintercept = 0, color = "darkgray", linetype = 2) +
    facet_wrap(~index, nrow = 1, scales = "free_x") +
    theme_bw() +
    labs(color = "", y = "") +
    theme(legend.position = "none")
  
  
  beta_tests_cat <-
    read_tsv("results/tables/permanova_host_vars_cat.tsv", col_types = cols()) %>% 
    filter(dataset %in% ds) %>% 
    filter(test_var == term) %>% 
    left_join(variables %>% select(test_var = var_id, var_name, variable_category = dataset, lididem_crc), by = join_by(test_var)) %>% 
    mutate(test_level = rename_var_levels(test_level)) %>% 
    mutate(var_level = paste(var_name, test_level)) %>% 
    rename(p = `Pr(>F)`)
  
  
  beta_assoc_plot <-
    beta_tests_cat %>% 
    mutate(dist = "Bray-Curtis") %>% 
    filter(!test_var %in% c("final_result", "detect_worthy_lesions")) %>% 
    filter(lididem_crc |
             str_detect(test_var, "Antibiotic") |
             var_name %in% c("WCRF", "Screening center", "National background", 
                             "Employment status", "Education", "Antibiotics use",
                             "Diverticulitis", "Hemorrhoids", "IBD"),
           test_level != "Missing") %>%
    ## Filter out levels with less than 10 observations
    filter(paste(term, test_level) %in% tmp_var_levels) %>% 
    group_by(test_var) %>% 
    mutate(m_R2 = mean(R2[ !test_level %in% c("Unknown", "Missing")])) %>% 
    ungroup() %>% 
    arrange(m_R2, R2) %>% 
    mutate(var_level = factor(var_level, levels = var_level)) %>% 
    ggplot(aes(x = R2, y = var_level, fill = -log10(p))) +
    geom_col(color = "black") +
    theme_bw() +
    facet_wrap(~dist) +
    scale_fill_gradient2(high = scales::muted("red"),
                         mid = "gray",
                         low = "white",
                         midpoint = 1) +
    labs(y = "")
  
  diversity_associations_plot <-
    ggpubr::ggarrange(alpha_assoc_plot + theme(text = element_text(size = 6)),
                      beta_assoc_plot + theme(text = element_text(size = 6), legend.position = "none"),
                      labels = c("d", "e"),
                      widths = c(1.2, 1))
  
  fig_1 <- ggpubr::ggarrange(prev_ab_plot,
                             diversity_associations_plot, 
                             nrow = 2, heights = c(1, 1.5))
  
  ggsave2("results/figures/microbiome_overview_figure.pdf", plot = fig_1, height = 150, width = 150, units = "mm")
  
}


plot_loc_associations <- function(ds = "metaphlan") {
  
  tmp_abund <- list("metaphlan" = mphlan_abundance,
                    "mags" = MAGs_abundance,
                    "keggs" = ko_hum)[[ds]]
  
  location_plot <-
    screening_data %>% 
    inner_join(sample_data %>% select(deltaker_id), by = join_by(deltaker_id)) %>% 
    select(deltaker_id, final_result) %>% 
    left_join(lesion_locs, by = "deltaker_id") %>% 
    filter(str_detect(lesion_diag, "^(5|6)")) %>% 
    mutate(lesion_diag = factor(lesion_diag, levels = c("5a. Non-advanced adenoma",
                                                        "5d. Advanced serrated lesion",
                                                        "5c. Advanced adenoma",
                                                        "6. Cancer"))) %>% 
    mutate(`Colonoscopy result` = factor(final_result, levels = c("5b. >= 3 Non-advanced adenomas",
                                                                  "5d. Advanced serrated",
                                                                  "5c. Advanced adenoma",
                                                                  "6. Cancer"))) %>% 
    (function(plot_obj) {
      plot_obj %>% 
        ggplot(aes(x = locX, y = -locY, color = lesion_diag)) +
        geom_point(data = plot_obj %>% filter(str_detect(lesion_diag, "^5a")), shape = 1, size = 0.6, alpha = 0.7) +
        geom_point(data = plot_obj %>% filter(str_detect(lesion_diag, "^5d")), shape = 1, size = 0.6, alpha = 0.7) +
        geom_point(data = plot_obj %>% filter(str_detect(lesion_diag, "^5c")), shape = 1, size = 0.6, alpha = 0.7) +
        geom_point(data = plot_obj %>% filter(str_detect(lesion_diag, "^6")), shape = 1, size = 0.6, alpha = 0.7) +
        facet_wrap(~`Colonoscopy result`, ncol = 1, labeller = label_wrap_gen(20)) +
        scale_color_manual(values = c(paletteer::paletteer_d(color_set)[c(5, 3:1, 9)])) +
        theme_bw() +
        coord_equal() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              text = element_text(size = 6),
              legend.position = "bottom") +
        labs(color = "Lesion diagnosis",
             x = "",
             y = "")
    })
  

  # Alpha diversity ---------------------------------------------------------

  alpha_diversity <- 
    tmp_abund %>% 
    pivot_longer(-sample_id) %>% 
    group_by(sample_id) %>% 
    summarize(Observed = sum(value > 0),
              Shannon = diversity(value),
              `Inverse Simpson` = diversity(value, index = "invsimpson"),
              .groups = "drop")
  
  loc_by_participant <-
    sample_data %>% 
    select(sample_id, deltaker_id) %>% 
    inner_join(screening_data %>% 
                 select(deltaker_id, most_serious_distancefromanus_processed, most_serious_segment, 
                        final_result, localization, distal_acn, proximal_acn,
                        sum_lesion_diameter_distal, sum_lesion_diameter_proximal, detect_worthy_lesions), 
               by = join_by(deltaker_id)) %>% 
    left_join(lesion_locs %>% 
                select(most_serious_segment = segment, 
                       segment_id, side_processed) %>% 
                distinct(), 
              by = join_by(most_serious_segment)) %>% 
    mutate(proximal_distal_lesion_summary = case_when(proximal_acn == 1 & distal_acn == 1 ~ "Both",
                                                      proximal_acn == 0 & distal_acn == 1 ~ "Distal",
                                                      proximal_acn == 1 & distal_acn == 0 ~ "Proximal",
                                                      proximal_acn == 0 & distal_acn == 0 ~ "None") %>% 
             factor(levels = c("None", "Proximal", "Distal", "Both")))
  
  a_div_dataset <-
    alpha_diversity %>% 
    pivot_longer(-sample_id, values_to = "Diversity", names_to = "Index") %>% 
    mutate(Index = factor(Index, levels = c("Observed", "Shannon", "Inverse Simpson"))) %>% 
    inner_join(loc_by_participant, join_by(sample_id))
  
  a_div_prox_dist_plot <-
    a_div_dataset %>% 
    filter(Index %in% "Observed") %>% 
    ggplot(aes(x = proximal_distal_lesion_summary, y = Diversity, fill = proximal_distal_lesion_summary)) +
    geom_jitter(width = 0.2, size = 0.1) +
    scale_fill_manual(values = c(paletteer::paletteer_d("awtools::bpalette")[c(10,7,2,4)])) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    facet_wrap(~"Presence of advanced colorectal lesions", nrow = 3, scales = "free_y") +
    labs(x = "Presence of advanced colorectal lesions",
         y = "Observed species") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 6))
  
  a_div_dist_plot <-
    a_div_dataset %>% 
    filter(Index %in% "Observed") %>% 
    filter(detect_worthy_lesions %in% "positive") %>%
    mutate(side_processed = factor(side_processed, labels = c("Distal localization", "Proximal localization"))) %>% 
    filter(!is.na(side_processed)) %>% 
    ggplot(aes(x = most_serious_distancefromanus_processed, y = Diversity)) +
    theme_bw() +
    geom_point(size = 0.3) +
    facet_wrap(~side_processed, ncol = 2) +
    labs(x = "Insertion depth from anal verge (cm)", 
         color = "Lesion localization",
         y = "Observed species") +
    theme(text = element_text(size = 6)) +
    stat_smooth(method = "lm")
  
  a_div_lesion_burden <-
    a_div_dataset %>% 
    filter(Index %in% "Observed") %>% 
    filter(detect_worthy_lesions %in% "positive") %>%
    pivot_longer(c(sum_lesion_diameter_distal, sum_lesion_diameter_proximal), 
                 names_to = "side", values_to = "sum_lesion_diameter") %>% 
    mutate(side = case_when(str_detect(side, "proximal") ~ "Proximal lesions",
                            str_detect(side, "distal") ~ "Distal lesions")) %>% 
    filter(sum_lesion_diameter != 0) %>% 
    ggplot(aes(x = sum_lesion_diameter, y = Diversity)) +
    geom_point(size = 0.4) +
    scale_x_log10() +
    facet_wrap(~side, ncol = 2) +
    stat_smooth(method = "lm") +
    theme_bw() +
    theme(text = element_text(size = 6)) +
    labs(x = "Sum of lesion diameters (mm)",
         y = "Observed species")
  
  tmp_dist <- 
    tmp_abund %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "sample_id") %>% 
    as.matrix() %>% 
    vegan::vegdist()
  
  tmp_pcoa_tab <-
    tmp_dist %>% 
    dist_to_PCoA(group_var = tmp_abund %>% 
                   select(sample_id) %>% 
                   left_join(loc_by_participant %>% select(sample_id, proximal_distal_lesion_summary), 
                             by = join_by(sample_id)) %>% 
                   pull(proximal_distal_lesion_summary))
  
  plot_pcoa_lesion_loc <- 
    tmp_pcoa_tab %>% 
    mutate(group = factor(group, levels = rev(levels(factor(group))))) %>% 
    plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2", type = "dot", with_bars = FALSE) +
    scale_color_manual(values = c(paletteer::paletteer_d("awtools::bpalette")[c(10,7,2,4)])) +
    scale_fill_manual(values = c(paletteer::paletteer_d("awtools::bpalette")[c(10,7,2,4)])) +
    theme(text = element_text(size = 6)) 
  
  
  
  localization_da_res <- 
    read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_localization_models.tsv", col_types = cols()) %>% 
    mutate(outcome = case_when(outcome == "lesion_dist_distal" ~ "Distance to distal lesion",
                               outcome == "lesion_dist_proximal" ~ "Distance to proximal lesion",
                               outcome == "lesion_loc_overall" ~ "Presence of advanced lesion",
                               outcome == "lesion_loc_size_distal" ~ "Sum of distal lesion diameters",
                               outcome == "lesion_loc_size_proximal" ~ "Sum of proximal lesion diameters") %>% 
             factor(levels = c("Presence of advanced lesion",
                               "Distance to distal lesion", 
                               "Distance to proximal lesion", 
                               "Sum of distal lesion diameters", 
                               "Sum of proximal lesion diameters"))) %>% 
    mutate(value = case_when(str_detect(value, "most_serious_distance") ~ "Insertion depth from anal verge",
                             str_detect(value, "distal_acn") ~ "Presence of distal advanced colorectal lesion",
                             str_detect(value, "proximal_acn") ~ "Presence of proximal advanced colorectal lesion",
                             str_detect(value, "diameter_distal") ~ "Sum of distal lesion diameters",
                             str_detect(value, "diameter_proximal") ~ "Sum of proximal lesion diameters")) %>% 
    mutate(significance = case_when(qval < 0.05 ~ "FDR significant",
                                    pval < 0.05 ~ "nominally significant",
                                    TRUE ~ "Not sigificant"))
  
  lesion_loc_da_plot <- 
    localization_da_res %>% 
    filter(str_detect(model, "wcrf.*antibiotics")) %>% 
    group_by(feature) %>% 
    filter(any(qval < 0.05)) %>% 
    ungroup() %>% 
    (function(x) {
      tmp <- 
        x %>% 
        filter(value %in% "Presence of distal advanced colorectal lesion") %>% 
        mutate(feature = fct_reorder(feature, desc(coef))) %>% 
        pull(feature) %>% 
        levels()
      
      x %>% 
        mutate(feature = factor(feature, levels = tmp))
    }) %>% 
    ggplot(aes(x = coef, y = feature, color = value)) +
    geom_vline(xintercept = 0, linetype = 2, color = "dark gray") +
    geom_point(aes(shape = significance), position = position_dodge2(width = 0.3)) +
    facet_wrap(~outcome, nrow = 1, labeller = label_wrap_gen(15)) +
    scale_y_discrete(labels = function(x) x %>% enframe(value = "sgb") %>% 
                       left_join(mp4_taxonomy, by = "sgb") %>% pull(species)) +
    scale_color_manual(values = paletteer::paletteer_d("fishualize::Balistapus_undulatus")[c(5:1,6)]) +
    scale_shape_manual(values = c(19,20,1)) +
    theme_bw() +
    labs(y = "",
         x = "log2FC",
         color = "") +
    theme(text = element_text(size = 6),
          axis.text.y = element_text(face = "italic", size = 6))
  
  a_div_loc_plot <- ggpubr::ggarrange(a_div_prox_dist_plot, 
                                      a_div_dist_plot, 
                                      a_div_lesion_burden, 
                                      ncol = 1, labels = c("b", "c", "d"), align = "h")
  
  localization_plot <-
    ggpubr::ggarrange(ggarrange(location_plot, NULL, heights = c(1, 0.1), ncol = 1, legend = "none"), 
                      ggpubr::ggarrange(a_div_loc_plot, 
                                        plot_pcoa_lesion_loc, 
                                        ncol = 1, labels = c("", "e"), 
                                        heights = c(2.5,1.5), legend = "none"), 
                      lesion_loc_da_plot, 
                      ncol = 3, labels = c("a", "", "f"), 
                      legend = "none", widths = c(1, 1, 2.5))
  
  ggsave2("results/figures/localization_plot.pdf", plot = localization_plot, height = 150, width = 210, units = "mm")
  
  
  ## Supplementary alpha
  
  a_div_prox_dist_plot_suppl <-
    a_div_dataset %>% 
    filter(!Index %in% "Observed") %>% 
    ggplot(aes(x = proximal_distal_lesion_summary, y = Diversity, fill = proximal_distal_lesion_summary)) +
    geom_jitter(width = 0.2, size = 0.1) +
    scale_fill_manual(values = c(paletteer::paletteer_d("awtools::bpalette")[c(10,7,2,4)])) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    facet_wrap(~Index, nrow = 2, scales = "free_y") +
    labs(x = "Presence of advanced colorectal lesions") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 6))
  
  a_div_dist_plot_suppl <-
    a_div_dataset %>% 
    filter(!Index %in% "Observed") %>% 
    filter(detect_worthy_lesions %in% "positive") %>%
    mutate(side_processed = factor(side_processed, labels = c("Distal localization", "Proximal localization"))) %>% 
    ggplot(aes(x = most_serious_distancefromanus_processed, y = Diversity)) +
    theme_bw() +
    geom_point(size = 0.3) +
    facet_wrap(~Index+side_processed, scales = "free_y", ncol = 2) +
    labs(x = "Insertion depth from anal verge (cm)", color = "Lesion localization") +
    theme(text = element_text(size = 6)) +
    stat_smooth(method = "lm")
  
  a_div_lesion_burden_suppl <-
    a_div_dataset %>% 
    filter(!Index %in% "Observed") %>% 
    filter(detect_worthy_lesions %in% "positive") %>%
    pivot_longer(c(sum_lesion_diameter_distal, sum_lesion_diameter_proximal), 
                 names_to = "side", values_to = "sum_lesion_diameter") %>% 
    mutate(side = case_when(str_detect(side, "proximal") ~ "Proximal lesions",
                            str_detect(side, "distal") ~ "Distal lesions")) %>% 
    filter(sum_lesion_diameter != 0) %>% 
    ggplot(aes(x = sum_lesion_diameter, y = Diversity)) +
    geom_point(size = 0.4) +
    scale_x_log10() +
    facet_wrap(~Index+side, scales = "free_y", ncol = 2) +
    stat_smooth(method = "lm") +
    theme_bw() +
    theme(text = element_text(size = 6)) +
    labs(x = "Sum of lesion diameters (mm)")
  
  a_div_loc_plot_suppl <- ggpubr::ggarrange(a_div_prox_dist_plot_suppl, 
                                            a_div_dist_plot_suppl, 
                                            a_div_lesion_burden_suppl, 
                                            ncol = 3, labels = "auto", align = "h")
  
  ggsave2("results/figures/localization_alpha_suppl.pdf", plot = a_div_loc_plot_suppl, height = 80,  width = 180, units = "mm")
  
}

plot_main_da_loc_assoc <- function(ds = "metaphlan") {
  
  if (!"ds" %in% ls()) {
    ds <- "metaphlan"
  }
  if (ds %in% "ko_hum") {
    ds_da <- "keggs"
  } else {
    ds_da <- ds
  }
  
  localization_da_res <- 
    read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_localization_models.tsv", col_types = cols()) %>% 
    mutate(outcome = case_when(outcome == "lesion_dist_distal" ~ "Distance to distal lesion",
                               outcome == "lesion_dist_proximal" ~ "Distance to proximal lesion",
                               outcome == "lesion_loc_overall" ~ "Presence of advanced lesion",
                               outcome == "lesion_loc_size_distal" ~ "Sum of distal lesion diameters",
                               outcome == "lesion_loc_size_proximal" ~ "Sum of proximal lesion diameters") %>% 
             factor(levels = c("Presence of advanced lesion",
                               "Distance to distal lesion", 
                               "Distance to proximal lesion", 
                               "Sum of distal lesion diameters", 
                               "Sum of proximal lesion diameters"))) %>% 
    mutate(value = case_when(str_detect(value, "most_serious_distance") ~ "Insertion depth from anal verge",
                             str_detect(value, "distal_acn") ~ "Presence of distal advanced colorectal lesion",
                             str_detect(value, "proximal_acn") ~ "Presence of proximal advanced colorectal lesion",
                             str_detect(value, "diameter_distal") ~ "Sum of distal lesion diameters",
                             str_detect(value, "diameter_proximal") ~ "Sum of proximal lesion diameters")) %>% 
    mutate(significance = case_when(qval < 0.05 ~ "FDR significant",
                                    pval < 0.05 ~ "nominally significant",
                                    TRUE ~ "Not sigificant")) %>% 
    filter(str_detect(model, "wcrf.*anti")) ## Standard model
  
  
  diff_abund_table <- 
    read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_main_models.tsv", col_types = cols()) %>% 
    mutate(value = str_replace(value, "positive", "CRC-related findings"))  %>%
    mutate(metadata = factor(metadata, 
                             levels = c("detect_worthy_lesions", "final_result_cat_neg"), 
                             labels = enframe(c("detect_worthy_lesions", "final_result_cat_neg"), value = "var_id") %>% 
                               left_join(variables, by = "var_id") %>% pull(var_name))) %>% 
    mutate(value = rename_var_levels(value)) %>% 
    filter(dataset %in% ds_da) %>% 
    filter(str_detect(model, "wcrf.*anti")) %>% ## Standard model
    group_by(feature) %>% 
    filter(any(qval < 0.05)) %>% 
    ungroup() %>% 
    (function(x) {
      x %>% 
        mutate(lab = factor(feature, 
                            levels = x %>% 
                              filter(value == "CRC-related findings") %>% 
                              arrange(desc(coef)) %>% 
                              pull(feature)))
    })
  
  main_da_loc_plot <- 
    localization_da_res %>% 
    inner_join(diff_abund_table %>% select(feature, lab) %>% distinct(), by = "feature") %>% 
    ggplot(aes(x = coef, y = lab, color = value)) +
    geom_vline(xintercept = 0, linetype = 2, color = "dark gray") +
    geom_point(aes(shape = significance), position = position_dodge2(width = 0.3)) +
    facet_wrap(~outcome, nrow = 1, labeller = label_wrap_gen(15)) +
    scale_y_discrete(labels = function(x) x %>% enframe(value = "sgb") %>% 
                       left_join(mp4_taxonomy, by = "sgb") %>% pull(species)) +
    scale_color_manual(values = paletteer::paletteer_d("fishualize::Balistapus_undulatus")[c(5:1,6)]) +
    scale_shape_manual(values = c(19,20,1)) +
    theme_bw() +
    labs(y = "",
         x = "log2FC",
         color = "") +
    theme(text = element_text(size = 6),
          axis.text.y = element_text(face = "italic", size = 6))
  
  ggsave(file = "results/figures/diff_abund_features_by_lesion_location_size.pdf", plot = main_da_loc_plot, height = 110, width = 180, units = "mm")
}

plot_microbial_outcome_associations <- function(ds = "metaphlan") {
  
  tmp_var_levels <-
    define_prev_var_levels(var_name_and_id = "id")
  
  if (!"ds" %in% ls()) {
    ds <- "metaphlan"
  }
  if (ds %in% "ko_hum") {
    ds_da <- "keggs"
  } else {
    ds_da <- ds
  }
  
  tmp_abund <- list("metaphlan" = mphlan_abundance,
                    "mags" = MAGs_abundance,
                    "keggs" = ko_hum)[[ds]]
  
  alpha_adj_outcome_tests_multinom <-
    read_tsv(file.path("results","tables","alpha_adj_outcome_tests.tsv"), col_types = cols()) %>% 
    filter(dataset %in% ds) %>% 
    filter(variable %in% "diversity") %>% 
    filter(model %in% "standard_model") %>% 
    filter(subset_cat %in% "all") %>% 
    filter(test_var != "final_result") %>% 
    mutate(index = factor(index, levels = c("observed", "shannon", "invsimpson"))) %>% 
    left_join(variables %>% select(var_name, test_var = var_id), by = "test_var") %>% 
    mutate(lvl_lab = case_when(y.level %in% "1" ~ "CRC-related findings",
                               TRUE ~ y.level)) %>% 
    mutate(lvl_lab = fct_relevel(lvl_lab, "CRC-related findings",
                                 levels(meta_dat$final_result_cat_neg)[-1]))
  
  alpha_outcome_clin_rel <-
    alpha_adj_outcome_tests_multinom %>% 
    filter(test_var %in% c("detect_worthy_lesions")) %>%
    ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = lvl_lab, color = lvl_lab)) +
    geom_vline(xintercept = 1, linetype = 2, color = "gray") +
    geom_pointrange(shape = 20) +
    facet_wrap(~index, ncol = 3) +
    theme_bw() +
    theme() +
    # scale_color_manual(values = paletteer::paletteer_d(color_set)[c(6)]) +
    scale_color_manual(values = color_assignments) +
    labs(x = "OR", 
         y = "",
         color = "") +
    labs(x = "") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 6))
  
  alpha_outcome_fin_res <-
    alpha_adj_outcome_tests_multinom %>% 
    filter(!test_var %in% c("detect_worthy_lesions")) %>%
    ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = lvl_lab, color = lvl_lab)) +
    geom_vline(xintercept = 1, linetype = 2, color = "gray") +
    geom_pointrange(shape = 20) +
    facet_wrap(~index, ncol = 3) +
    theme_bw() +
    theme() +
    scale_color_manual(values = color_assignments) +
    # scale_color_manual(values = paletteer::paletteer_d(color_set)[c(5:1)]) +
    labs(x = "OR", 
         y = "",
         color = "") +
    theme(strip.text = element_blank(),
          text = element_text(size = 6))
  
  lims <- rbind(ggplot_build(alpha_outcome_clin_rel)$layout$panel_scales_x[[1]]$range$range, 
                ggplot_build(alpha_outcome_fin_res)$layout$panel_scales_x[[1]]$range$range) %>% 
    (function(x) {
      c(min(x[,1]), max(x[,2]))
    })
  
  plot_alpha_outcome_assoc <-
    ggpubr::ggarrange(alpha_outcome_clin_rel + scale_x_continuous(limits = lims),
                      alpha_outcome_fin_res + scale_x_continuous(limits = lims), 
                      ncol = 1, heights = c(0.35, 1),
                      align = "v",
                      legend = "none")
  
  tmp_dist <- 
    tmp_abund %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "sample_id") %>% 
    as.matrix() %>% 
    vegan::vegdist()
  
  tmp_pcoa_tab <-
    tmp_dist %>% 
    dist_to_PCoA(group_var = tmp_abund %>% 
                   select(sample_id) %>% 
                   left_join(sample_data %>% select(sample_id, deltaker_id), by = "sample_id") %>% 
                   left_join(screening_data %>% select(deltaker_id, final_result_cat_neg), by = "deltaker_id") %>% 
                   pull(final_result_cat_neg))
  
  plot_pcoa_fin_res <- 
    tmp_pcoa_tab %>% 
    mutate(group = factor(group, levels = levels(meta_dat$final_result_cat_neg))) %>%
    inner_join(sample_data %>% select(sample_id, deltaker_id), by = "sample_id") %>% 
    left_join(screening_data %>% select(deltaker_id, final_result_cat_neg), by = "deltaker_id") %>%
    plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2", type = "dot", with_bars = TRUE) %>% 
    lapply(function(x) {
      if (is.null(x)) x
      else {
        x + 
          scale_color_manual(values = color_assignments) +
          scale_fill_manual(values = color_assignments) +
          theme(legend.position = "none",
                plot.margin = unit(rep(0,4), "lines"),
                text = element_text(size = 6))
      }
    }) %>% 
    ggpubr::ggarrange(plotlist = ., align = "hv", 
                      widths = c(1,0.25), heights = c(0.2, 1)) +
    coord_fixed()
  
  dmm_group <-
    load_dmm(synthetic_data = TRUE) %>% 
    mutate(gr = rename_var_levels(gr))
  
  tmp_pcoa_tab_dmm <-
    tmp_dist %>% 
    dist_to_PCoA(group_var = tmp_abund %>% 
                   select(sample_id) %>% 
                   left_join(dmm_group %>% select(sample_id, gr), by = "sample_id") %>%
                   pull(gr))
  
  plot_pcoa_dmms <- 
    tmp_pcoa_tab_dmm %>% 
    plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2", type = "dot", with_bars = TRUE) %>% 
    lapply(function(x) {
      if (is.null(x)) x
      else {
        x + 
          scale_color_manual(values = color_assignments) +
          scale_fill_manual(values = color_assignments) +
          theme(legend.position = "none",
                plot.margin = unit(rep(0,4), "lines"),
                text = element_text(size = 6))
      }
    }) %>% 
    ggpubr::ggarrange(plotlist = ., align = "hv", 
                      widths = c(1,0.25), heights = c(0.2, 1)) +
    coord_fixed()
  
  
  dmm_mnom_outcome_assoc_one_v_rest <-
    read_tsv("results/tables/dmm/dmm_mnom_outcome_assoc_one_v_rest.tsv", col_types = cols()) %>% 
    filter(term %in% "grtest_level") %>% 
    mutate(y.level = str_replace(y.level, "^1$", "CRC-related findings")) %>% 
    mutate(dmm_test_level = rename_var_levels(dmm_test_lvl)) %>% 
    arrange(dmm_test_lvl) %>% 
    mutate(dmm_test_level = factor(dmm_test_level,
                                   levels = unique(dmm_test_level), 
                                   labels = paste0(unique(dmm_test_level)," versus others"))) %>% 
    mutate(lvl_lab = fct_relevel(y.level, "CRC-related findings",
                                 levels(meta_dat$final_result_cat_neg)[-1]))
  
  dmm_outcome_clin_rel <-
    dmm_mnom_outcome_assoc_one_v_rest %>% 
    filter(outcome %in% c("detect_worthy_lesions")) %>%
    ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = lvl_lab, color = lvl_lab)) +
    geom_vline(xintercept = 1, linetype = 2, color = "gray") +
    geom_pointrange(shape = 20) +
    facet_wrap(~dmm_test_level, ncol = 4) +
    theme_bw() +
    scale_color_manual(values = color_assignments) +
    labs(x = "OR", 
         y = "",
         color = "") +
    labs(x = "") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 6)) +
    coord_fixed() +
    scale_x_continuous(trans = "log10")
  
  dmm_outcome_fin_res <-
    dmm_mnom_outcome_assoc_one_v_rest %>% 
    filter(!outcome %in% c("detect_worthy_lesions")) %>%
    ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = lvl_lab, color = lvl_lab)) +
    geom_vline(xintercept = 1, linetype = 2, color = "gray") +
    geom_pointrange(shape = 20) +
    facet_wrap(~dmm_test_level, ncol = 4) +
    theme_bw() +
    scale_color_manual(values = color_assignments) +
    labs(x = "OR", 
         y = "",
         color = "") +
    theme(strip.text = element_blank(),
          text = element_text(size = 6)) +
    scale_x_continuous(trans = "log10")
  
  dmm_lims <- 
    rbind(ggplot_build(dmm_outcome_clin_rel)$layout$panel_scales_x[[1]]$range$range, 
          ggplot_build(dmm_outcome_fin_res)$layout$panel_scales_x[[1]]$range$range) %>% 
    (function(x) {
      10^c((min(x[,1])), max(x[,2]))
    })
  
  plot_dmm_outcome_assoc <-
    ggpubr::ggarrange(dmm_outcome_clin_rel + coord_cartesian(xlim = dmm_lims),
                      dmm_outcome_fin_res + coord_cartesian(xlim = dmm_lims), 
                      ncol = 1, heights = c(0.35, 1),
                      align = "v",
                      legend = "none")
  
  
  diversity_plot <-
    ggarrange(plot_alpha_outcome_assoc, 
              plot_pcoa_fin_res, 
              nrow = 1,
              labels = c("a", "b"))
  
  dmm_plots <- ggarrange(plot_pcoa_dmms,
                         plot_dmm_outcome_assoc, 
                         nrow = 1,
                         labels = c("c", "d"))
  
  diversity_dmm_outcome_plot <-
    ggarrange(diversity_plot, dmm_plots, nrow = 2)
  
  ggsave2(paste0("results/figures/diversity_association_",ds,".pdf"), plot = diversity_dmm_outcome_plot, height = 160, width = 220, units = "mm")
}


plot_pairwise_beta_div_final_res <- function() {
  
  pairwise_fin_res_dist_comp <-
    read_tsv("results/tables/beta_fin_res_pairwise.tsv", col_types = cols())
  
  pairwise_beta_plot <-
    pairwise_fin_res_dist_comp %>% 
    filter(term %in% "final_result_cat_neg") %>% 
    mutate(across(c(test_level,ref_level), .fns = function(x) factor(x, levels = levels(meta_dat$final_result_cat_neg)))) %>% 
    filter(as.integer(ref_level) < as.integer(test_level)) %>% 
    filter(adjustment %in% "standard_model") %>% 
    ggplot(aes(x = ref_level, y = test_level, fill = F, label = ifelse(`Pr(>F)` < 0.05, "*", ""))) +
    geom_tile() +
    geom_text() +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          text = element_text(size = 6)) +
    labs(x = "Reference level",
         y = "Test level",
         fill = "Pseudo-F")
  
  ggsave2(filename = "results/figures/pairwise_beta.pdf", plot = pairwise_beta_plot, height = 80, width = 100, units = "mm")
}

plot_pairwise_diff_abund <- function() {
  
  pairwise_diff_abund_res <-
    read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_pairwise_models.tsv", col_types = cols()) %>% 
    mutate(value = str_remove(value, "^(0|1)"))
  
  tmp <- 
    pairwise_diff_abund_res %>% 
    group_by(feature) %>% 
    filter(any(qval < 0.05)) %>% 
    ungroup() %>% 
    mutate(across(.cols = c(value, ref), function(x) factor(x) %>% fct_relevel(levels(meta_dat$final_result_cat_neg)))) %>% 
    (function(x) {
      x %>% 
        mutate(feature = factor(feature, 
                                levels = x %>% 
                                  filter(str_detect(ref, "Neg"), str_detect(value, "Canc")) %>% 
                                  arrange(desc(coef)) %>% 
                                  pull(feature)))  
    })
    
  
  tmp_abund <-
    mphlan_abundance %>% 
    pivot_longer(-sample_id, values_to = "abundance", names_to = "sgb") %>% 
    inner_join(tmp %>% select(sgb = feature) %>% distinct(), by = join_by(sgb)) %>% 
    left_join(sample_data %>% select(sample_id, deltaker_id), by = "sample_id") %>% 
    left_join(meta_dat %>% select(deltaker_id, final_result_cat_neg), by = "deltaker_id") %>% 
    group_by(sgb) %>% 
    mutate(scaled_abundance = scale(abundance)) %>% 
    group_by(final_result_cat_neg, sgb) %>% 
    summarize(mean_abundance = mean(scaled_abundance), .groups = "drop") %>% 
    mutate(sgb = factor(sgb, levels = levels(tmp$feature)))
  
  tmp2 <-
    tmp %>% 
    filter(qval < 0.05) %>% 
    filter(as.integer(ref) < as.integer(value))
  
  pairwise_da_plot <- 
    tmp_abund %>% 
    ggplot(aes(x = final_result_cat_neg, y = sgb)) +
    geom_tile(aes(fill = mean_abundance)) +
    scale_y_discrete(labels = function(x) x %>% enframe(value = "sgb") %>% left_join(mp4_taxonomy, by = "sgb") %>% pull(species)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = "Diagnostic group",
         y = "Species",
         fill = "Mean of scaled abundance") +
    scale_fill_gradient2() +
    geom_segment(aes(x = ref, xend = value, 
                     y = as.integer(feature)+as.integer(ref)*0.1-0.2, 
                     yend = as.integer(feature)+as.integer(ref)*0.1-0.2),
                 data = tmp2) +
    theme(axis.text.y = element_text(face = "italic"))
  
  ggsave("results/figures/pairwise_da.pdf", pairwise_da_plot, height = 8, width = 6.5)
    
}

summarize_dmm <- function(k = 4, dataset = "crcbiome") {
  
  posterior_mean_diff <- 
    read_tsv("results/tables/dmm/posterior_mean_diff.tsv", col_types = cols())
  
  if (dataset == "crcbiome") {
    tmp_abund <- mphlan_abundance
  }
  
  lab_func <- function(labs_in, names = "sgb", lab_out = "species", df = mp4_taxonomy) {
    labs_in %>% 
      enframe(value = names) %>% 
      left_join(df, by = names) %>% 
      pull(any_of(lab_out))
  }
  
  tax_contrib_dirich_comp <-
    posterior_mean_diff %>% 
    select(sgb, all_of(paste0("V",seq(k)))) %>% 
    slice_head(n = 30) %>% 
    mutate(sgb = factor(sgb, levels = rev(sgb))) %>% 
    pivot_longer(-sgb, values_to = "prob", names_to = "group") %>% 
    mutate(group = rename_var_levels(group))
  
  tax_contrib_plot <-
    tax_contrib_dirich_comp %>% 
    ggplot(aes(x = group, y = sgb, fill = prob)) +
    scale_y_discrete(labels = lab_func) +
    geom_tile() +
    scale_fill_viridis_c(trans = "log10") +
    theme_bw() +
    labs(fill = "Contribution",
         y = "",
         x = "Component") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ## DMM associations
  mnom_group_host_var_assoc <-
    read_tsv("data/dmm/dmm_host_assoc_lrt.tsv", col_types = cols()) %>% 
    filter(!var_id %in% "bowel_disorder_merged_2")
  
  ## Plot
  dmm_host_assoc_plot <- 
    mnom_group_host_var_assoc %>% 
    filter(!is.na(Test)) %>% 
    mutate(model_complexity = factor(model_complexity, 
                                     levels = c("host_assoc_model", "standard_model"))) %>% 
    arrange(desc(`Pr(Chi)`)) %>% 
    left_join(variables %>% select(var_id, var_name), by = "var_id") %>% 
    mutate(var_name = factor(var_name, levels = unique(var_name[ model_complexity %in% "host_assoc_model"]))) %>% 
    filter(!is.na(var_name), !str_detect(var_id, "(detect_worthy_lesions|final_result)"), model_complexity == "host_assoc_model") %>% 
    ggplot(aes(x = -log10(`Pr(Chi)`), y = var_name, fill = model_complexity)) +
    geom_col(position = "dodge") +
    geom_vline(xintercept = -log10(0.05), color = "red", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 5),
          legend.position = "none") +
    scale_fill_manual(values = c(paletteer::paletteer_d(color_set)[c(2:3)])) +
    labs(y = "",
         fill = "",
         x = "-log10(LRT p-value)") 
  
  
  dmm_plot <-
    ggpubr::ggarrange(tax_contrib_plot + theme(text = element_text(size = 8)),
                      dmm_host_assoc_plot, labels = c("a", "b"), nrow = 1, widths = c(1, 1.2))
  
  dmm_plot %>% 
    ggsave(filename = "results/figures/dmm_plots.pdf", plot = ., height = 100, width = 190, units = "mm")
  
}



summarize_diff_abund <- function(ds = "metaphlan") {
  
  if (!"ds" %in% ls()) {
    ds <- "metaphlan"
  }
  if (ds %in% "ko_hum") {
    ds_da <- "keggs"
  } else {
    ds_da <- ds
  }
  
  tmp_abund <- list("metaphlan" = mphlan_abundance,
                    "mags" = MAGs_abundance,
                    "keggs" = ko_hum)[[ds_da]]
  
  tmp_var_levels <-
    define_prev_var_levels()
  
  diff_abund_table <- 
    read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_main_models.tsv", col_types = cols()) %>%
    mutate(value = str_replace(value, "positive", "CRC-related findings"))  %>%
    mutate(metadata = factor(metadata, 
                             levels = c("detect_worthy_lesions", "final_result_cat_neg"), 
                             labels = enframe(c("detect_worthy_lesions", "final_result_cat_neg"), value = "var_id") %>% 
                               left_join(variables, by = "var_id") %>% pull(var_name))) %>% 
    mutate(value = rename_var_levels(value)) %>% 
    filter(dataset %in% ds_da) %>% 
    filter(str_detect(model, "wcrf.*anti")) ## Standard model
  
  da_host_assoc <- 
    read_tsv("data/diff_abund/Maaslin2/meta_vars/maaslin2_meta_vars.tsv", col_types = cols()) %>% 
    mutate(dataset = case_when(dataset %in% "mphlan" ~ "metaphlan",
                               dataset %in% "ko_hum" ~ "keggs",
                               dataset %in% "MAGs" ~ "mags"))     %>% 
    left_join(variables %>% 
                select(variable = var_id, var_name), 
              by = "variable")
  
  tmp_da_dataset <-
    diff_abund_table %>% 
    filter(!str_detect(value, "4.")) %>%
    group_by(feature) %>% 
    mutate(any_sign = any(qval[ (metadata %in% "Outcome") | 
                                  (metadata %in% "CRC-related findings")] < 0.05)) %>% 
    filter(any_sign) %>% 
    ungroup() 
  
  if (ds_da %in% "metaphlan") {
    lab_tab <- mp4_taxonomy %>% 
      select(feature = sgb, lab = species)
  } else if (ds_da %in% "keggs") {
    lab_tab <-
      ko_annotations %>% 
      select(feature = KO_ID, lab = KO_name)
  } else if (ds_da %in% "mags") {
    lab_tab <- 
      MAG_taxonomy %>% 
      select(feature = MAG_id, lab = species)
  }
  
  ## Add labels
  tmp_da_dataset2 <- 
    tmp_da_dataset %>% 
    left_join(lab_tab, by = "feature") %>% 
    mutate(feature = factor(feature, levels = tmp_da_dataset %>% 
                              filter(metadata %in% "CRC-related findings") %>% 
                              arrange(desc(coef)) %>% 
                              pull(feature))) %>% 
    mutate(lab = case_when(nchar(lab) > 40 ~ paste(str_extract(lab, ".{40}"), "..."),
                           TRUE ~ lab)) %>% 
    mutate(metadata = factor(metadata, 
                             levels = c("CRC-related findings", "Outcome")))
  
  tmp_da_dataset3 <-
    tmp_da_dataset2 %>%
    left_join(read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_unadjusted_models.tsv", col_types = cols()) %>% 
                mutate(value = str_replace(value, "positive", "CRC-related findings")) %>% 
                filter(dataset %in% ds_da) %>%
                select(feature, value, unadj_coef = coef, unadj_qval = qval), 
              by = join_by(feature, value))
  

  main_outcome_da_plot <-
    tmp_da_dataset2 %>% 
    ggplot(aes(x = coef, y = feature, color = value)) +
    geom_vline(xintercept = 0, linetype = 2, color = "dark gray") +
    geom_pointrange(aes(xmin = coef-stderr, xmax = coef+stderr, alpha = qval < 0.05), 
                    position = position_dodge(width = 0.65), 
                    size = 0.2,
                    shape = 1) +
    geom_point(aes(x = unadj_coef), 
               data = tmp_da_dataset3 %>% filter(!metadata %in% "Outcome", unadj_qval < 0.05), 
               shape = 20, alpha = 1, size = 0.4, color = "red") +
    facet_wrap(~metadata, labeller = label_wrap_gen(width = 15)) +
    scale_y_discrete(labels = function(x) x %>% enframe(value = "feature") %>% left_join(tmp_da_dataset2 %>% select(feature, lab) %>% distinct(), by = "feature") %>% pull(lab)) +
    theme_bw() +
    scale_alpha_manual(values = c(0.3, 1)) +
    scale_color_manual(values = color_assignments) +
    labs(x = "log2FC", y = "")
  
  
  ## Prevalence/abundance
  ab_prev_summary <-
    tmp_abund %>% 
    pivot_longer(-sample_id, names_to = "feature", values_to = "abundance") %>% 
    inner_join(tmp_da_dataset2 %>% select(feature) %>% unique(), by = "feature") %>% 
    group_by(feature) %>% 
    summarize(mean_abundance = ifelse(ds == "metaphlan", mean(abundance/100), mean(abundance)),
              prevalence = sum(abundance > 0)/n()) %>% 
    mutate(feature = factor(feature, levels = levels(tmp_da_dataset2$feature))) 
  
  prev_sum_plot <-
    ab_prev_summary %>% 
    ggplot(aes(x = prevalence, y = feature)) +
    geom_col() +
    scale_x_continuous(breaks = c(0,1), limits = c(0,1)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 6)) +
    labs(x = "", y = "") +
    facet_wrap(~" ")
  
  ab_sum_plot <-
    ab_prev_summary %>% 
    ggplot(aes(x = mean_abundance, y = feature)) +
    geom_point() +
    scale_x_continuous(trans = "log10")  +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 6)) +
    labs(x = "", y = "") +
    facet_wrap(~" ")

  
  prev_abund_summary_plot <-
    ggpubr::ggarrange(prev_sum_plot,
                      ab_sum_plot,
                      nrow = 1, legend = "none",
                      align = "h")
  
  ## Host var associations
  tmp_outcome_host_da <- 
    da_host_assoc %>% 
    mutate(x = "Participant characteristics") %>% 
    filter(dataset %in% ds_da) %>% 
    filter(feature %in% tmp_da_dataset2$feature) %>% 
    filter(!str_detect(variable, "(detect_worthy_lesions|final_result)")) %>% 
    mutate(feature = factor(feature, levels = levels(tmp_da_dataset2$feature))) %>% 
    mutate(sign_level = case_when(qval < 0.05 ~ sprintf('\u2217'),
                                  TRUE ~ "")) %>% 
    left_join(variables %>% select(variable = var_id, var_name, variable_category = dataset, lididem_crc), by = join_by(variable, var_name)) %>% 
    mutate(test_level = str_remove(value, "^[ab]_")) %>% 
    filter(paste(var_name, test_level) %in% c(define_prev_var_levels(), define_prev_var_levels(name_update = FALSE))) %>% 
    group_by(var_name) %>% 
    filter(any(qval[!str_detect(test_level, "(Unknown|Missing)") & 
                      !str_detect(var_name, "(FIT|Initial|round)")] < 0.05)) %>% 
    ungroup() 

  tmp_outcome_host_da <-
    tmp_outcome_host_da %>% 
    mutate(var_name = factor(var_name, levels = (tmp_outcome_host_da %>% 
                                                   group_by(feature, var_name) %>% 
                                                   summarize(mean_coef = mean(coef), .groups = "drop") %>% 
                                                   order_by_cluster_variable("var_name", "feature", "mean_coef")))) %>% 
    arrange(var_name) %>% 
    mutate(test_level = rename_var_levels(test_level)) %>% 
    mutate(var_level = paste(var_name, test_level)) %>% 
    mutate(var_level = factor(var_level, levels = unique(var_level)))
  
  
  
  da_host_assoc_plot <-
    tmp_outcome_host_da %>% 
    ggplot(aes(x = var_level, y = feature, fill = coef, label = sign_level)) +
    geom_tile() +
    geom_text(size = 2) +
    scale_fill_gradient2() +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1), 
          ) +
    facet_wrap(~x) +
    labs(x = "", y = "", fill = "log2FC") +
    guides(fill = guide_colourbar(theme = theme(
      legend.key.width  = unit(.5, "lines"),
      legend.key.height = unit(2.5, "lines")
    )))
  
  diff_abund_plot <-
    ggpubr::ggarrange(main_outcome_da_plot + 
                        theme(legend.position = "none", 
                              text = element_text(size = 6),
                              axis.text.y = element_text(size = 6, face = "italic")),
                      prev_sum_plot + theme(legend.position = "none"),
                      ab_sum_plot + theme(legend.position = "none"),
                      da_host_assoc_plot + 
                        theme(text = element_text(size = 6),
                              axis.text.x = element_text(size = 6, angle = 45),
                              axis.ticks.y = element_blank()),
                      align = "h", widths = c(1,0.2, 0.2,1.5), nrow = 1,
                      labels = c("a", "b", "", "c"))
  
  if (ds %in% "metaphlan") {
    ## To enable unicode asterisks
    cairo_pdf(paste0("results/figures/diff_abund.pdf"), height = 120*0.03937, width = 250*0.03937)
    diff_abund_plot
    dev.off()
  } else {
    ## To enable unicode asterisks
    cairo_pdf(paste0("results/figures/diff_abund_",ds,".pdf"), height = 120*0.03937, width = 250*0.03937)
    diff_abund_plot
    dev.off()
  }
}

plot_strat_da <- function(ds_da = "metaphlan") {
  
  diff_abund_table <- 
    read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_main_models.tsv", col_types = cols()) %>% 
    mutate(value = str_replace(value, "positive", "CRC-related findings"))  %>%
    mutate(metadata = factor(metadata, 
                             levels = c("detect_worthy_lesions", "final_result_cat_neg"), 
                             labels = enframe(c("detect_worthy_lesions", "final_result_cat_neg"), value = "var_id") %>% 
                               left_join(variables, by = "var_id") %>% pull(var_name))) %>% 
    mutate(value = rename_var_levels(value)) %>% 
    filter(dataset %in% ds_da) %>% 
    filter(str_detect(model, "wcrf.*anti")) ## Standard model
  
  da_res_strat <- 
    read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_strat_models_fix.tsv", col_types = cols())
    
  sensitivity_da_plot <-
    da_res_strat %>% 
    filter((strat_var %in% c("antibiotics_reg_quest_comb")) |
             strat_var %in% "PPI_antacids_reg_quest_comb" |
             (strat_var %in% "colo_any")) %>%
    left_join(variables %>% select(strat_var = var_id, var_name) %>% bind_rows(tibble(strat_var = "colo_any", var_name = "Non-neoplastic findings")), by = join_by(strat_var)) %>% 
    mutate(significance = case_when(qval < 0.05 ~ "FDR significant",
                                    pval < 0.05 ~ "Nominally significant",
                                    TRUE ~ "Not significant") %>% 
             factor(levels = c("Not significant", "Nominally significant", "FDR significant"))) %>% 
    mutate(strat_val_lab = paste0(strat_val, " (n=", N, ")")) %>% 
    mutate(value = ifelse(value == "positive", "CRC-related findings", value)) %>% 
    inner_join(diff_abund_table %>% 
                 group_by(feature) %>% 
                 filter(any(qval < 0.05)) %>%
                 ungroup() %>% 
                 (function(x) {
                   x %>% 
                     mutate(lab = factor(feature, 
                                         levels = x %>% 
                                           filter(value == "CRC-related findings") %>% 
                                           arrange(desc(coef)) %>% 
                                           pull(feature)))
                 }) %>% 
                 filter(qval < 0.05) %>%
                 select(feature, lab, value) %>% 
                 distinct(), by = join_by(feature, value)) %>% 
    mutate(value = factor(value, levels = c("CRC-related findings", levels(meta_dat$final_result_cat_neg)[-1]))) %>% 
    (function(plot_obj) {
      plot_obj %>% 
        ggplot(aes(x = coef, y = lab, color = value)) +
        geom_vline(xintercept = 0, linetype = 2, color = "dark gray") +
        geom_point(aes(shape = significance), 
                   position = position_dodge(width = 0.65)) +
        scale_shape_manual(values = c(1,20,19)) +
        facet_wrap(~var_name + strat_val_lab, labeller = label_wrap_gen(width = 15), nrow = 1) +
        scale_y_discrete(labels = function(x) x %>% enframe(value = "sgb") %>% left_join(mp4_taxonomy, by = "sgb") %>% pull(species)) +
        theme_bw() +
        # scale_alpha_manual(values = c(0.3, 1)) +
        scale_color_manual(values = color_assignments) +
        labs(x = "log2FC", y = "", color = "") +
        theme(text = element_text(size = 6),
              axis.text.y = element_text(size = 6, face = "italic"))
    })
  
  ggsave2("results/figures/diff_abund_sensitivity.pdf", sensitivity_da_plot, height = 100, width = 180, units = "mm")
}


plot_ecoli_lit_rep_da <- function() {
  
  # E coli ------------------------------------------------------------------
  ecoli_pres_abs_models <-
    read_tsv("results/tables/ecoli_pks/ecoli_pks_detection_outcomes.tsv", col_types = cols()) %>% 
    mutate(test_level = factor(test_level, levels = c("CRC-related findings", levels(meta_dat$final_result_cat_neg)[-1]))) %>% 
    mutate(pks_strat = factor(pks_strat, labels = c("E coli detection overall", "E coli detection stratified by pks status")))
  
  ecoli_outcome_crc_rel <-
    ecoli_pres_abs_models %>% 
    filter(outcome_var %in% c("detect_worthy_lesions")) %>%
    ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = test_level, color = test_level, shape = pks_lab)) +
    geom_vline(xintercept = 1, linetype = 2, color = "gray") +
    geom_pointrange(position = position_dodge(width = 0.3)) +
    facet_wrap(~pks_strat) +
    theme_bw() +
    theme() +
    scale_color_manual(values = color_assignments) +
    labs(x = "OR", 
         y = "",
         color = "") +
    labs(x = "") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 6))
  
  ecoli_outcome_fin_res <-
    ecoli_pres_abs_models %>% 
    filter(outcome_var %in% c("final_result_cat_neg")) %>%
    ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = test_level, color = test_level, shape = pks_lab)) +
    geom_vline(xintercept = 1, linetype = 2, color = "gray") +
    geom_pointrange(position = position_dodge(width = 0.3)) +
    facet_wrap(~pks_strat) +
    theme_bw() +
    theme() +
    scale_color_manual(values = color_assignments) +
    labs(x = "OR", 
         y = "",
         color = "") +
    labs(x = "") +
    theme(strip.text = element_blank(),
          text = element_text(size = 6))
  
  lims <- rbind(ggplot_build(ecoli_outcome_crc_rel)$layout$panel_scales_x[[1]]$range$range, 
                ggplot_build(ecoli_outcome_fin_res)$layout$panel_scales_x[[1]]$range$range) %>% 
    (function(x) {
      c(min(x[,1]), max(x[,2]))
    })
  
  plot_ecoli_outcome_assoc <-
    ggpubr::ggarrange(ecoli_outcome_crc_rel + scale_x_continuous(limits = lims, trans = "log10"),
                      ecoli_outcome_fin_res + scale_x_continuous(limits = lims, trans = "log10"), 
                      ncol = 1, heights = c(0.35, 1),
                      align = "v",
                      legend = "none")
  
  # Lit rep associations ----------------------------------------------------
  
  crc_species <- 
    read_tsv("data/misc/piccino.csv", col_types = cols()) %>% 
    mutate(reported_enrichment = ifelse(Effect_size_adj > 0, "CRC", "Control")) %>% 
    arrange(desc(Effect_size_adj))
  
  da_main_low_prev <- 
    read_tsv("data/diff_abund/Maaslin2/main_outcome/maaslin2_summary_low_prev_models.tsv", col_types = cols()) %>% 
    mutate(metadata = factor(metadata, 
                             levels = c("detect_worthy_lesions", "final_result_cat_neg"), 
                             labels = enframe(c("detect_worthy_lesions", "final_result_cat_neg"), value = "var_id") %>% 
                               left_join(variables, by = "var_id") %>% pull(var_name))) %>% 
    filter(dataset %in% "metaphlan") %>% 
    select(-c(dataset, strat)) %>% 
    left_join(crc_species %>% select(feature = sgb_alt, reported_enrichment) %>% distinct(), by = "feature") %>% 
    mutate(reported_enrichment = case_when(is.na(reported_enrichment) ~ "no enrichment",
                                           TRUE ~ reported_enrichment)) %>% 
    mutate(prevalence_cat = ifelse(N.not.zero/N > 0.1, "high prevalence", "low prevalence")) %>% 
    mutate(feature_lab = factor(feature, 
                                levels = (crc_species %>% filter(sgb_alt %in% mp4_taxonomy$sgb) %>% select(sgb_alt) %>% distinct() %>%  pull(sgb_alt)),
                                labels = (crc_species %>% 
                                            filter(sgb_alt %in% mp4_taxonomy$sgb) %>% 
                                            distinct(sgb_alt, .keep_all = TRUE) %>% 
                                            mutate(lab = paste(Microbial_feature, sgb_alt)) %>% 
                                            pull(lab)))) %>%
    mutate(value = fct_relevel(value, "positive") %>% fct_relabel(.fun = function(x) ifelse(x == "positive", "CRC-related findings", x))) %>% 
    filter(str_detect(metadata, "(CRC-related|Outcome)")) %>% 
    mutate(significance = case_when(qval < 0.05 ~ "FDR significant",
                                    pval < 0.05 ~ "Nominally significant",
                                    TRUE ~ "Not significant") %>% 
             factor(levels = c("Not significant", "Nominally significant", "FDR significant")))
  
  lit_assoc_es_plot <-
    da_main_low_prev %>% 
    filter(reported_enrichment != "no enrichment") %>% 
    group_by(feature) %>% 
    filter(any(pval[ str_detect(value, "^Cancer")] < 0.05)) %>% 
    ungroup() %>% 
    (function(x) {
      n_CRC <- x %>% filter(reported_enrichment == "Cancer") %>% pull(feature) %>% unique() %>% length()
      n_control <- x %>% filter(reported_enrichment == "Control") %>% pull(feature) %>% unique() %>% length()
      new_lvls <- tibble(feature_lab = paste("x", sapply(seq(n_CRC-n_control), function(y) paste(sample(letters, 3, replace = TRUE), collapse = ""))),
                         reported_enrichment = "Control")
      x %>% 
        bind_rows(new_lvls) %>% 
        mutate(feature_lab = factor(feature_lab, levels = c(levels(x$feature_lab), new_lvls$feature_lab)))
      
    }) %>% 
    
    ggplot(aes(x = coef, y = feature_lab, color = value)) +
    geom_vline(xintercept = 0, linetype = 2, color = "dark gray") +
    geom_point(aes(shape = significance), 
               position = position_dodge(width = 0.65)) +
    scale_shape_manual(values = c(1,20,19)) +
    facet_wrap(~reported_enrichment, scales = "free_y", labeller = label_wrap_gen(width = 15), ncol = 1) +
    theme_bw() +
    scale_color_manual(values = color_assignments) +
    scale_x_continuous(minor_breaks = NULL) +
    labs(x = "log2FC", y = "", color = "") +
    theme(text = element_text(size = 6),
          axis.text.y = element_text(size = 6, face = "italic"))
  
  
  ecoli_lit_rep_fig <- 
    ggarrange(plot_ecoli_outcome_assoc,
              lit_assoc_es_plot, ncol = 1, heights = c(1, 1.6),
              legend = "none",
              labels = "auto")
  
  ggsave2("results/figures/ecoli_lit_rep_associations.pdf", 
          plot = ecoli_lit_rep_fig, height = 200, width = 80, units = "mm")
  

  # Additional lit rep plots ------------------------------------------------

  lit_assoc_volcano_plot <-
    da_main_low_prev %>%
    ggplot(aes(x = coef, y = -log10(pval), color = reported_enrichment, shape = pval < 0.05)) +
    geom_point(data = da_main_low_prev %>% filter(reported_enrichment %in% "no enrichment"), color = "gray", alpha = 0.4) +
    geom_point(data = da_main_low_prev %>% filter(!reported_enrichment %in% "no enrichment")) +
    scale_shape_manual(values = c(1,20)) +
    facet_wrap(~prevalence_cat + value, nrow = 2) +
    theme_bw() +
    labs(x = "log2FC") +
    theme(text = element_text(size = 6))

  ggsave2("results/figures/lit_rep_associations_volcano.pdf", plot = lit_assoc_volcano_plot, height = 100, width = 220, units = "mm")

  lit_assoc_full_es_plot <-
    da_main_low_prev %>%
    filter(reported_enrichment != "no enrichment") %>%
    ggplot(aes(x = coef, y = feature_lab, color = reported_enrichment)) +
    geom_vline(xintercept = 0, color = "gray", linetype = 2) +
    geom_point(aes(shape = pval < 0.05)) +
    scale_shape_manual(values = c(1,20)) +
    facet_wrap(~value, nrow = 1) +
    theme_bw() +
    labs(x = "log2FC") +
    # scale_y_discrete(labels = function(x) x %>% enframe(value = "sgb") %>% left_join(crc_species, by = "sgb") %>% pull(species)) +
    theme(text = element_text(size = 6),
          axis.text.y = element_text(size = 4, face = "italic"))

  ggsave2("results/figures/lit_rep_associations_plot_all_species.pdf", plot = lit_assoc_full_es_plot, height = 220, width = 180, units = "mm")
}

algorithm_and_dataset_comparison <- function(primary_model = "mphl", secondary_alg = "log_reg", secondary_add_var = "") {
  
  alg_and_dataset_comp <-
    read_tsv("data/models/evaluations/alg_comp_models/primary_assessments.tsv", col_types = cols()) %>% 
    rename(model = primary_model,
           primary_dataset = dataset) %>% 
    left_join(make_model_summary(.), by = "model") %>%
    filter(!is.na(dataset) & !data_type %in% "mag-based")
  
  (plot_alg_and_dataset_comp <-
      alg_and_dataset_comp %>%
      filter(predictor %in% "primary_prob",
             algorithm != "neural net") %>% 
      group_by(model, data_type, dataset, model_name, algorithm, dataset_type) %>%
      summarize_metric(aggregation = "by_rep") %>%
      ungroup() %>%
      ggplot(aes(x = auroc_m_mean, xmin = auroc_m_mean-auroc_sem, xmax = auroc_m_mean+auroc_sem, y = dataset, color = algorithm)) +
      geom_vline(xintercept = 0.5, linetype = 2, color = "red") +
      geom_pointrange(position = position_dodge(width = 0.65), size = 0.1) +
      paletteer::scale_color_paletteer_d(color_set) +
      geom_vline(xintercept = 0.5, linetype = 2, color = "gray") +
      scale_x_continuous(limits = c(0.45, 0.7)) +
      ggstats::geom_stripped_rows(linewidth = 0) +
      theme_bw() +
      labs(x = "Mean CV AUROC",
           color = "Algorithm",
           y = "Training dataset") +
      theme(text = element_text(size = 6))
  )
  
  ggsave2("results/figures/ml_alg_dataset_comp.pdf", 
          plot = plot_alg_and_dataset_comp, 
          height = 80,
          width = 110,
          units = "mm")
  
}


summarize_ml_res <- function(ds = "metaphlan") {
  
  tmp_models <- read_tsv("data/models.tsv", col_types = cols())
  
  if (!"ds" %in% ls()) {
    ds <- "metaphlan"
  }
  if (ds %in% "metaphlan") {
    tmp_dat <- "mphl"
  }
  
  process_main_model <- function(mod_res, part = "training") {
    
    ## Fix variable name
    mod_res <-
      mod_res %>% 
      rename(any_of(c(primary_dataset = "dataset")))
    
    mod_res %>% 
      filter(primary_dataset %in% "mphl") %>% 
      filter((secondary_model == "log_reg_mphl_fit" & predictor %in% c("primary_prob", "FIT_value", "secondary_prob")) |
               (secondary_model == "log_reg_mphl_fit_demo_wcrf" & predictor %in% c("secondary_prob", "only_add_var"))) %>% 
      mutate(mod = case_when(predictor %in% "FIT_value" ~ "FIT alone",
                             predictor %in% "primary_prob" ~ "microbial model alone",
                             predictor %in% "secondary_prob" & sec_add_var %in% "fit" ~ "microbial model + FIT",
                             predictor %in% "secondary_prob" & sec_add_var %in% "fit_demo_wcrf" ~ "microbial model + FIT, sex, age, and WCRF score",
                             predictor %in% "only_add_var" ~ "FIT, sex, age, and WCRF score") %>% 
               factor(levels = c("FIT alone", "microbial model alone", "microbial model + FIT",
                                 "microbial model + FIT, sex, age, and WCRF score",
                                 "FIT, sex, age, and WCRF score")))
  }
  
  outcome_strat_models <-
    read_tsv("data/models/evaluations/main_models/assessments.tsv", col_types = cols()) %>% 
    process_main_model() %>% 
    mutate(selection = factor(selection, levels = c("all", "neg-v-naa", "neg-v-as", "neg-v-aa", "neg-v-crc"))) %>% 
    group_by(mod, selection) %>% 
    summarize_model_stats(stats = c("auroc", "specificity", "sensitivity"), output_wide = FALSE) %>% 
    mutate(plot_measure = case_when(measure == "auroc" ~ "AUROC",
                                    measure %in% "sensitivity" ~ "Sensitivity",
                                    measure %in% "specificity" ~ "Specificity")) %>% 
    mutate(xint = ifelse(plot_measure == "AUROC", 0.5, NA))
  
  ## Standard model - evaluation of identification of distal lesions - should probably not include those with only proximal lesions
  distal_models <-
    read_tsv("data/models/evaluations/main_models/distal_assessments.tsv", col_types = cols()) %>% 
    process_main_model() %>% 
    group_by(mod) %>% 
    summarize_model_stats(stats = c("auroc", "specificity", "sensitivity"), output_wide = FALSE) 
  
  proximal_models <-
    read_tsv("data/models/evaluations/main_models/proximal_assessments.tsv", col_types = cols()) %>% 
    process_main_model() %>% 
    group_by(mod) %>% 
    summarize_model_stats(stats = c("auroc", "specificity", "sensitivity"), output_wide = FALSE) 
  
  ## Variable strat predictions
  strat_models <- 
    read_tsv("data/models/evaluations/stratified_predictions/assessments.tsv", col_types = cols()) %>% 
    process_main_model() %>% 
    group_by(mod, strat_var, strat_val, n_in_group) %>% 
    summarize_model_stats(stats = c("auroc"), output_wide = FALSE) %>% 
    mutate(strat_val = rename_var_levels(strat_val)) %>%
    filter(!(strat_var %in% "FIT_value_cat" & mod %in% "FIT alone")) 
  strat_models %>% filter(strat_var == "dmm")
  
  main_mod_outcome_enum <-
    read_tsv("data/models/evaluations/main_models/outcome_enumerations_tmp.tsv", col_types = cols()) %>% 
    process_main_model() %>% 
    group_by(rep, fold, mod, enum_var, threshold, fraction_def_var) %>% 
    summarize(`Classified positive` = sum(n[ prediction])/sum(n),
              `Classified negative` = sum(n[ !prediction])/sum(n), 
              .groups = "drop") %>% 
    group_by(mod, enum_var, threshold, fraction_def_var) %>% 
    summarize_model_stats(stats = "Classified positive", output_wide = FALSE) %>% 
    mutate(perc_reduced = (1-threshold)*100) %>% 
    mutate(fraction_def_var = rename_var_levels(fraction_def_var)) %>% 
    mutate(enum_var = case_when(enum_var == "negative" ~ "Not CRC-related",
                                enum_var == "positive" ~ "CRC-related",
                                TRUE ~ enum_var))
  
  
  limit_points <- 
    outcome_strat_models %>% 
    filter(selection %in% "all",
           mod %in% "microbial model alone") %>% 
    mutate(mean = case_when(measure %in% "auroc" ~ 0.5,
                            measure %in% "sensitivity" ~ 0.5,
                            measure %in% "specificity" ~ 0.1)) %>% 
    bind_rows(outcome_strat_models %>% 
                filter(selection %in% "all",
                       mod %in% "microbial model alone") %>% 
                mutate(mean = case_when(measure %in% "auroc" ~ 0.7,
                                        measure %in% "sensitivity" ~ 1,
                                        measure %in% "specificity" ~ 0.65)))
  
  ## Test set eval
  fix_test_input <- function(df) {
    df %>%
      filter(str_detect(primary_model, tmp_dat)) %>% 
      mutate(predictor = case_when(predictor %in% "secondary_pred" & no_m ~ "only_add_var",
                                   TRUE ~ predictor) %>% 
               str_replace("pred", "prob")) %>% 
      rename(sec_add_var = add_var) %>% 
      mutate(dataset = tmp_dat) %>% 
      mutate(secondary_model = paste0("log_reg_", dataset, "_", sec_add_var)) %>% 
      process_main_model()
  }
  
  test_aurocs <- 
    read_tsv("data/models/final_model/evaluations/AUROC.tsv", col_types = cols()) %>% fix_test_input()
  final_roc_obj <-
    read_tsv("data/models/final_model/evaluations/roc.tsv", col_types = cols()) %>% fix_test_input()
  final_sens_spec <-
    read_tsv("data/models/final_model/evaluations/sens_spec.tsv", col_types = cols()) %>% 
    mutate(predictor = str_replace(predictor, "_pred", "_prob")) %>% 
    mutate(primary_dataset = str_remove(primary_model, "rf_")) %>% 
    rename(sec_add_var = add_var) %>% 
    mutate(secondary_model = paste0("log_reg_", str_remove(primary_model, "rf_"), "_", sec_add_var)) %>% 
    process_main_model() #%>%
  
  
  ## Plots
  
  (plot_training_classification_main_models <-
      outcome_strat_models %>%
      filter(selection %in% c("all")) %>% 
      ggplot(aes(x = mean,
                 xmin = mean-sd,
                 xmax = mean+sd,
                 color = mod,
                 y = mod)) +
      geom_vline(aes(xintercept = xint), 
                 data = outcome_strat_models %>% 
                   filter(!is.na(xint)), 
                 linetype = 2, color = "darkgrey") +
      geom_pointrange(shape = 20, size = 0.25) +
      geom_point(alpha = 0, data = limit_points) +
      scale_color_manual(values = color_assignments) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(color = "",
           y = "",
           x = "Mean CV model performance") +
      theme(text = element_text(size = 8)) +
      facet_wrap(~plot_measure, ncol = 3, scales = "free_x"))
  
  (roc_final_model <-
      final_roc_obj %>% 
      filter(var_cat %in% "roc") %>% 
      ggplot(aes(x = specificities, y = sensitivities, group = mod, color = mod)) +
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed") +
      geom_path() +
      scale_color_manual(values = color_assignments) +
      coord_fixed() +
      scale_x_reverse() +
      theme_bw() +
      theme(text = element_text(size = 8)) +
      labs(color = "")
  )
  
  
  (plot_stratified_classification_main <-
      strat_models %>% 
      filter(strat_var %in% c("dmm", "kjonn", "senter", "WCRF", "age_cat", "FIT_value_cat")) %>%
      mutate(strat_var = factor(strat_var,
                                levels = c("dmm", "kjonn", "age_cat", "senter", "WCRF", "FIT_value_cat"),
                                labels = c("Community state", "Sex", "Age", "Screening center", "WCRF", "FIT value"))) %>%
      mutate(strat_val = fct_relevel(strat_val, "low", "mid", "high")) %>% 
      arrange(mod, strat_var, strat_val) %>% 
      mutate(lab = sprintf("%s \n(n = %3s)", strat_val, n_in_group)) %>% 
      mutate(lab = factor(lab, levels = unique(lab))) %>% 
      ggplot(aes(x = mean, xmin = mean-sd, xmax = mean+sd, y = lab, color = mod)) +
      geom_vline(xintercept = 0.5, linetype = 2, color = "darkgrey") +
      geom_pointrange(position = position_dodge(width = 0.65), shape = 20, size = 0.25) +
      ggstats::geom_stripped_rows(linewidth = 0) +
      scale_color_manual(values = color_assignments) +
      scale_x_continuous(limits = c(0.4, NA)) +
      theme_bw() +
      labs(color = "",
           x = "Mean CV AUROC",
           y = "") +
      theme(text = element_text(size = 8)) +
      facet_wrap(~strat_var, ncol = 2, scales = "free_y"))
  
  
  (plot_training_classification_outcome_comparison <-
      outcome_strat_models %>%
      filter(selection != "all") %>% 
      filter(measure %in% c("auroc", "sensitivity")) %>% 
      ggplot(aes(x = mean,
                 xmin = mean-sd,
                 xmax = mean+sd,
                 color = mod,
                 y = selection)) +
      geom_vline(aes(xintercept = xint), 
                 data = outcome_strat_models %>% 
                   filter(!is.na(xint)), 
                 linetype = 2, color = "darkgrey") +
      geom_pointrange(position = position_dodge(width = 0.5), shape = 20, size = 0.25) +
      ggstats::geom_stripped_rows(linewidth = 0) +
      scale_color_manual(values = color_assignments) +
      theme_bw() +
      labs(color = "",
           y = "",
           x = "Mean CV model performance") +
      theme(text = element_text(size = 8)) +
      facet_wrap(~plot_measure, ncol = 3, scales = "free_x"))
  
  
  (plot_loc_spec_classification_comparison <-
      distal_models %>%
      mutate(selection = "Negative vs advanced distal lesions") %>% 
      bind_rows(proximal_models %>%
                  mutate(selection = "Negative vs advanced proximal lesions")) %>% 
      filter(measure %in% "auroc") %>% 
      ggplot(aes(x = mean,
                 xmin = mean-sd,
                 xmax = mean+sd,
                 color = mod,
                 y = selection)) +
      geom_vline(xintercept = 0.5, linetype = 2, color = "darkgrey") +
      geom_pointrange(position = position_dodge(width = 0.5), shape = 20, size = 0.25) +
      scale_x_continuous(limits = c(0.45, NA)) +
      scale_y_discrete(labels = scales::wrap_format(20)) +
      ggstats::geom_stripped_rows(linewidth = 0) +
      scale_color_manual(values = color_assignments) +
      theme_bw() +
      theme(text = element_text(size = 8)) +
      labs(color = "",
           y = "", 
           x = "Mean CV AUROC") +
      facet_wrap(~""))
  
  
  
  ml_plot <-
    ggarrange(
      ggarrange(plot_training_classification_main_models,
                roc_final_model,
                ncol = 2,
                legend = "none",
                labels = c("a", "b")),
      ggarrange(plot_stratified_classification_main + theme(legend.position = "none"),
                ggarrange(plot_training_classification_outcome_comparison,
                          plot_loc_spec_classification_comparison, 
                          ncol = 1, 
                          legend = "none",
                          labels = c("d", "e")),
                ncol = 2, 
                labels = "c",
                widths = c(1.5, 1)),
      ncol = 1,
      heights = c(1, 1.5))
  
  ggsave2("results/figures/ml_summary.pdf", plot = ml_plot, height = 180, width = 180, units = "mm")
  
  
  ## Suppl figures
  
  (main_mod_outcome_enum_plot <-
      main_mod_outcome_enum %>% 
      mutate(fraction_def_var = rename_var_levels(fraction_def_var)) %>% 
      ggplot(aes(y = mean, x = perc_reduced, group = enum_var, color = enum_var)) +
      geom_line() +
      geom_point(size = 0.3) +
      facet_wrap(~mod+fraction_def_var, ncol = 2, labeller = label_wrap_gen(width = 25)) +
      labs(x = "Percent reduction in number of colonoscopy referrals",
           y = "Mean CV fraction of group classified as positive",
           color = "Colonoscopy-based outcomes") +
      theme_bw() +
      theme(text = element_text(size = 6)) +
      scale_color_manual(values = color_assignments))
  
  ggsave2("results/figures/enumeration_by_outcome.pdf", plot = main_mod_outcome_enum_plot, height = 150, width = 100, units = "mm")
  
  (plot_stratified_classification_suppl <-
      strat_models %>% 
      mutate(strat_var = case_when(strat_var == "colo_hemorrhoids" ~ "Hemorrhoids",
                                   strat_var == "colo_diverticulitis" ~ "Diverticulitis",
                                   TRUE ~ strat_var)) %>% 
      filter(strat_var %in% c("Cardiovascular disease", "Chronic obstructive pulmonary disease", 
                              "Diabetes", "antibiotics_use", 
                              "Smoking", "Utdanning",
                              "Hemorrhoids", "Diverticulitis")) %>%
      mutate(strat_var = factor(strat_var,
                                levels = c("Smoking", 
                                           "Utdanning",
                                           "antibiotics_use",
                                           "Hemorrhoids",
                                           "Cardiovascular disease", 
                                           "Chronic obstructive pulmonary disease", 
                                           "Diabetes",
                                           "Diverticulitis")) %>% 
               fct_recode(`Antibiotics use` = "antibiotics_use",
                          Education = "Utdanning")) %>% 
      # labels = c("Enterotype", "Sex", "Screening center", "WCRF", "Antibiotics use", "FIT value"))) %>%
      mutate(strat_val = fct_relevel(strat_val, "Primary school", "High school", "University/college", "No drug dispensed", "Drug dispensed")) %>% 
      arrange(mod, strat_var, strat_val) %>% 
      mutate(lab = sprintf("%s \n(n = %3s)", strat_val, n_in_group)) %>% 
      mutate(lab = factor(lab, levels = unique(lab))) %>% 
      ggplot(aes(x = mean, xmin = mean-sd, xmax = mean+sd, y = lab, color = mod)) +
      geom_vline(xintercept = 0.5, linetype = 2, color = "darkgrey") +
      geom_pointrange(position = position_dodge(width = 0.65), shape = 20, size = 0.25) +
      ggstats::geom_stripped_rows(linewidth = 0) +
      scale_color_manual(values = color_assignments) +
      theme_bw() +
      labs(color = "",
           x = "Mean CV AUROC",
           y = "") +
      facet_wrap(~strat_var, nrow = 2, scales = "free_y", labeller = label_wrap_gen(width = 20)))
  
  
  n_features_frac_samples_forests_selection <-
    read_tsv("data/models/n_features_frac_samples/mphlan_rf_selection.tsv", col_types = cols())
  
  n_features_frac_samples_plot <-
    n_features_frac_samples_forests_selection %>%
    group_by(n_vars_requested, frac_samples) %>% 
    summarize(mean_cv_auroc = mean(auroc),
              n_samples = mean(training_samples), 
              .groups = "drop") %>% 
    group_by(frac_samples) %>% 
    mutate(mean_n_samples = paste0(100*unique(frac_samples), "% (", round(mean(n_samples)), ")"), 
           .groups = "drop") %>% 
    mutate(mean_n_samples = fct_reorder(mean_n_samples, n_samples)) %>% 
    ggplot(aes(x = factor(mean_n_samples), 
               y = factor(n_vars_requested),
               label = signif(mean_cv_auroc, 2))) +
    geom_tile(aes(fill = mean_cv_auroc)) +
    geom_text(size = 1.5) +
    theme_bw() +
    labs(x = "Fraction of samples (n) used for training",
         y = "Number of features used for prediction",
         fill = "Mean CV AUROC") +
    theme(text = element_text(size = 6),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          legend.key.size = unit(4, units = "mm")) +
    coord_fixed()
  
  
  ml_res_plot_suppl <-
    ggpubr::ggarrange(
      plot_stratified_classification_suppl + theme(text = element_text(size = 8)),
      ggarrange(NULL, n_features_frac_samples_plot, labels = c("", "b")),
      ncol = 1, common.legend = TRUE, legend = "left", labels = c("a", ""), heights = c(1.2,1))
  
  ggsave2("results/figures/ml_suppl_fig.pdf",
          plot = ml_res_plot_suppl, height = 150, width = 220, units = "mm")
  
  
}


