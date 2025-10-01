

preprocess_input_data <- function() {
  
  # load data ---------------------------------------------------------------
  
  sample_meta <- read_rds("data/metadata/data_by_sample_20240321.Rds")
  
  screening_data <- read_rds("data/input_preprocessed/screening/20240212_screening_proc.rds")%>% 
    mutate(detect_worthy_lesions = case_when(final_result_cat7 %in% c("Advanced adenoma", "Advanced serrated", "CRC", 
                                                                      "Non-advanced adenoma (>=3)") ~ "positive",
                                             final_result_cat7 %in% c("Negative", "Non-advanced serrated/other lesion",
                                                                      "Non-advanced adenoma (<3)") ~ "negative")) %>% 
    mutate(age_cat = cut(age_invitation, breaks = c(0,67, 100), labels = c("67 or under", "over 67"))) %>% 
    tibble()
  
  prior_samples <- read_tsv("data/input_preprocessed/screening/20221108_FIT_prior_samples.tsv")
  
  first_round_invitation <- 
    read_csv("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/RAW/screening/mb_screening_all_rounds_2021-12-08.csv", 
             col_types = cols()) %>% 
    filter(runde == 1)
  
  lifestyle_data <- read_rds("data/input_preprocessed/lifestyle/241122_lifestyle_proc.Rds") %>% 
    mutate(`bowel_disorder_merged_2` = case_when(Bowel_disorder_merged %in% c("IBS", "IBD", "Celiac disease", "Other") ~ "Yes",
                                                                  Bowel_disorder_merged %in% "No bowel disease" ~ "No",
                                                                  TRUE ~ Bowel_disorder_merged)) %>% 
    tibble()
  
  diet_data <- read_rds("data/input_preprocessed/diet/241122_diet_proc.rds") %>% 
    tibble()
  
  lesions_data <- read_csv("/tsd/p1068/data/durable/007-f_smei/001-trro/CRCbiome/RAW/screening/2021_10/mb_lesions_2021-10-15.csv")
  lesions_size_data <- read_csv("data/input_preprocessed/Lesions_count_final.csv", col_types = cols()) %>% 
    select(diameter, Polyp_vol, deltaker_id, lesionno)
  
  lesions_size_by_participant <- 
    read_csv("data/input_preprocessed/Lesions_totalvolume_by_patients.csv", col_types = cols()) %>% 
    select(deltaker_id, 
           sum_polyp_volume = Total_polyp_vol,
           sum_proximal_polyp_volume = TotVol_Proximal,
           sum_distal_polyp_volume = TotVol_Distal)
  
  # screening_unproc <- read_csv("/tsd/p1068/data/durable/001-trro/CRCbiome/RAW/screening/2021_10/mb_screening_2021-10-15.csv")
  
  ## WCRF - Already excluded those with low quality questionnaires. Also excluded stage IV cancers (n = 2)
  diet_wcrf <- read_csv2("data/input_preprocessed/diet/151222_wcrf_proc_CRCbiome.csv")
  
  ## energy-adjusted foodgroups
  diet_ea <- read_csv2("data/input_preprocessed/diet/241122_nutrients_inc_suppl_residuals.csv")
  
  ## Metaphlan3
  header_names <- read_tsv("data/input_preprocessed/metagenome/tax_species.tsv", n_max = 1, col_names = FALSE) %>%
    pivot_longer(everything()) %>% 
    pull(value)
  
  abundance_data <- read_tsv("data/input_preprocessed/metagenome/tax_species.tsv", skip = 1, col_names = FALSE) %>%
    {if (length(header_names) == ncol(.)) {
      set_names(., header_names)
    } else {
      set_names(., c("feature", header_names))
    }} %>% 
    rename(feature = 1)
  
  ## Pivot
  abundance_data <- abundance_data %>% 
    pivot_longer(-feature) %>% 
    pivot_wider(names_from = feature, values_from = value) %>% 
    rename(sample_id = name)
  
  rm(header_names)
  
  read.delim("data/input_preprocessed/metagenome/tax_species_n.tsv") %>% 
    tibble() %>% 
    mutate(taxonomy_id = str_extract(clade_taxid, "(?<=/|)[:digit:]*$")) %>% 
    write_tsv("data/input_processed/metaphlan_species_taxonomy.tsv")
  
  ## Metaphlan 4
  mp4 <- read_tsv("data/input_preprocessed/profiles.tsv", skip = 1)
  
  mp4_sgbs <-
    mp4 %>% 
    filter(str_detect(clade_name, "t__")) %>% 
    mutate(clade_name = str_remove(clade_name, "^.*t__")) %>% 
    pivot_longer(-clade_name, names_to = "sample_id", values_to = "abundance") %>% 
    group_by(sample_id) %>% 
    mutate(abundance = abundance/sum(abundance)*100) %>% 
    ungroup() %>% 
    pivot_wider(names_from = clade_name, values_from = abundance)
  
  mp4_species <-
    mp4 %>% 
    filter(!str_detect(clade_name, "t__"),
           str_detect(clade_name, "s__")) %>% 
    mutate(clade_name = str_remove(clade_name, "\\|t__.*"),
           clade_name = str_remove(clade_name, ".*s__")) %>% 
    pivot_longer(-clade_name, names_to = "sample_id", values_to = "abundance") %>% 
    group_by(sample_id) %>% 
    mutate(abundance = abundance/sum(abundance)*100) %>% 
    ungroup() %>% 
    pivot_wider(names_from = clade_name, values_from = abundance)
  
  mp4_genera <-
    mp4 %>% 
    filter(!str_detect(clade_name, "s__"),
           str_detect(clade_name, "g__")) %>% 
    mutate(clade_name = str_remove(clade_name, "\\|s__.*"),
           clade_name = str_remove(clade_name, ".*g__")) %>% 
    pivot_longer(-clade_name, names_to = "sample_id", values_to = "abundance") %>% 
    group_by(sample_id) %>% 
    mutate(abundance = abundance/sum(abundance)*100) %>% 
    ungroup() %>% 
    pivot_wider(names_from = clade_name, values_from = abundance)
  
  mp4_sgb_names <- 
    mp4 %>% 
    filter(str_detect(clade_name, "t__")) %>% 
    select(clade_name) %>% 
    separate(clade_name, into = c("kingdom", "phylum", "clade", "order", "family", "genus", "species", "sgb"), 
             sep = "\\|", 
             remove = FALSE) %>% 
    mutate(across(-clade_name, .fns = function(x) x %>% str_remove("[:alpha:]__")))
  
  mp4_species_names <- 
    mp4 %>% 
    filter(!str_detect(clade_name, "t__"),
           str_detect(clade_name, "s__")) %>% 
    select(clade_name) %>% 
    separate(clade_name, into = c("kingdom", "phylum", "clade", "order", "family", "genus", "species"), 
             sep = "\\|", 
             remove = FALSE) %>% 
    mutate(across(-clade_name, .fns = function(x) x %>% str_remove("[:alpha:]__")))
  
  mp4_genera_names <- 
    mp4 %>% 
    filter(!str_detect(clade_name, "s__"),
           str_detect(clade_name, "g__")) %>% 
    select(clade_name) %>% 
    separate(clade_name, into = c("kingdom", "phylum", "clade", "order", "family", "genus"), 
             sep = "\\|", 
             remove = FALSE) %>% 
    mutate(across(-clade_name, .fns = function(x) x %>% str_remove("[:alpha:]__")))
  
  
  
  
  
  ## MAGs
  median_coverage_genomes <- read_tsv("data/input_preprocessed/metagenome/median_coverage_genomes.tsv")
  raw_counts_genomes <- read_tsv("data/input_preprocessed/metagenome/raw_counts_genomes.tsv")
  
  mags_taxonomy_db_metadata <- read_tsv("/ess/p1068/cluster/databases/atlas_db/GTDB_V05/metadata/genome_metadata.tsv") %>% 
    select(accession, ncbi_taxid, ncbi_taxonomy) %>% 
    rename(taxonomy_id = ncbi_taxid)
  
  mags_taxa_id <- read_tsv("data/input_preprocessed/metagenome/MAG_taxonomy/gtdbtk.bac120.summary.tsv") %>% 
    mutate(db = "bac120",
           closest_placement_radius = as.double(closest_placement_radius), 
           closest_placement_ani = as.double(closest_placement_ani), 
           closest_placement_af = as.double(closest_placement_af)) %>% 
    bind_rows(read_tsv("data/input_preprocessed/metagenome/MAG_taxonomy/gtdbtk.ar122.summary.tsv") %>% 
                mutate(db = "ar122")) %>% 
    rename(MAG_id = user_genome,
           accession = fastani_reference) %>% 
    mutate(species = str_extract(classification, "(?<=s__).*$"),
           genus = str_extract(classification, "(?<=g__).*(?=;s__)"),
           family = str_extract(classification, "(?<=f__).*(?=;g__)"),
           order = str_extract(classification, "(?<=o__).*(?=;f__)"),
           clade = str_extract(classification, "(?<=c__).*(?=;o__)"),
           phylum = str_extract(classification, "(?<=p__).*(?=;c__)"),
           domain = str_extract(classification, "(?<=d__).*(?=;p__)")) %>% 
    left_join(mags_taxonomy_db_metadata)
  
  
  ## MAG eggnog - full
  # MAG_eggnog <- read_tsv("data/input_preprocessed/metagenome/MAG_eggnog.tsv")
  ## MAG eggnog - COGs
  MAG_COG <- read_tsv("data/input_preprocessed/metagenome/MAG_COGs.tsv")
  ## MAG GO
  MAG_GO <- read_tsv("data/input_preprocessed/metagenome/MAG_GOs.tsv")
  ## MAG KO (KEGG orthology groups)
  MAG_KO <- read_tsv("data/input_preprocessed/metagenome/MAG_KOs.tsv")
  
  
  ## Bracken
  bracken_reads <- read_tsv("data/input_preprocessed/metagenome/bracken_reads.tsv")
  bracken_relabund <- read_tsv("data/input_preprocessed/metagenome/bracken_rel_abund.tsv")
  
  ## Viruses
  viral_abundance <- read_tsv("data/input_preprocessed/metagenome/coverage_table_min75.tsv")
  
  ## Pathways
  pathway_abundance <- read_tsv("data/input_preprocessed/metagenome/pathabundance.tsv")
  
  import_genedata <- 
    function(filename, obs_frac = 0.1, name_append = "_Abundance-RPKs") {
      read_tsv(filename) %>% 
        rename(ID = 1) %>% 
        rename_with(.fn = function(x) gsub(name_append, "", x)) %>% 
        select(c("ID", sample_meta %>% 
                   filter(grepl("S_", sample_id),
                          Total_Bases_QC_ATLAS >= 1e9,
                          Prøvetype %in% "Baseline") %>% 
                   pull(sample_id))) %>% 
        filter(apply(.[,-1]>0, 1, sum) > ((ncol(.)-1)*obs_frac)) %>% 
        filter(!ID %in% c("UNMAPPED", "UNGROUPED", "UNINTEGRATED")) %>% 
        pivot_longer(-ID, names_to = "sample_id") %>% 
        left_join(sample_meta %>% select(sample_id, reads_proc)) %>% 
        mutate(value = value/reads_proc*1e6) %>% 
        select(-reads_proc) %>% 
        pivot_wider(names_from = "ID", values_from = value) %>% 
        return()
    }
  
  ## humann data
  pathways <- import_genedata("data/input_preprocessed/metagenome/pathabundance.tsv", obs_frac = 0.1, name_append = "_Abundance")
  pfam <- import_genedata("data/input_preprocessed/metagenome/pfam.tsv")
  eggnog <- import_genedata("data/input_preprocessed/metagenome/eggnog.tsv")
  kegg <- import_genedata("data/input_preprocessed/metagenome/KEGG_orthogroups.tsv")
  go <- import_genedata("data/input_preprocessed/metagenome/go.tsv")
  l4ec <- import_genedata("data/input_preprocessed/metagenome/l4ec.tsv")
  ## Need different input function for genefamilies - too big for dplyr.
  # genefamilies <- import_genedata("data/input_preprocessed/metagenome/genefamilies.tsv")
  
  ## Annotation of ko's
  ko_annotation <- read_rds("data/input_preprocessed/metagenome/KO_function_annotation.Rds") %>% tibble()
  ko_annotation_h <- read_tsv("data/input_preprocessed/metagenome/humann_rerun/map_ko_name.txt", col_names = FALSE) %>% 
    dplyr::rename(ko_id = 1, ko_name = 2)
  
  # ko_annotation_h %>% 
  #   full_join(ko_annotation) %>% 
  #   mutate(a = is.na(KO_Name),
  #          b = is.na(KO_name)) %>% 
  #   count(a,b)
  
  ko_annotation %>% 
    write_tsv("data/input_processed/ko_annotations_w_brite.tsv")
  ko_annotation_h %>% 
    write_tsv("data/input_processed/humann_rerun/ko_annotations.tsv")
  
  pathways %>% write_rds("data/input_processed/pathways.Rds")
  pfam %>% write_rds("data/input_processed/pfam.Rds")
  eggnog %>% write_rds("data/input_processed/eggnog.Rds")
  kegg %>% write_rds("data/input_processed/kegg.Rds")
  go %>% write_rds("data/input_processed/go.Rds")
  l4ec %>% write_rds("data/input_processed/l4ec.Rds")
  
  
  ## humann data - rerun with mp4
  pathways <- import_genedata("data/input_preprocessed/metagenome/humann_rerun/pathabundance_no_tax.tsv", obs_frac = 0.1, name_append = "_Abundance")
  pfam <- import_genedata("data/input_preprocessed/metagenome/humann_rerun/pfam_no_tax.tsv")
  eggnog <- import_genedata("data/input_preprocessed/metagenome/humann_rerun/eggnog_no_tax.tsv")
  kegg <- import_genedata("data/input_preprocessed/metagenome/humann_rerun/ko_no_tax.tsv")
  go <- import_genedata("data/input_preprocessed/metagenome/humann_rerun/go_no_tax.tsv")
  l4ec <- import_genedata("data/input_preprocessed/metagenome/humann_rerun/level4ec_no_tax.tsv")
  ## Need different input function for genefamilies - too big for dplyr.
  # genefamilies <- import_genedata("data/input_preprocessed/metagenome/genefamilies.tsv")
  
  
  pathways %>% write_rds("data/input_processed/humann_rerun/pathways.Rds")
  pfam %>% write_rds("data/input_processed/humann_rerun/pfam.Rds")
  eggnog %>% write_rds("data/input_processed/humann_rerun/eggnog.Rds")
  kegg %>% write_rds("data/input_processed/humann_rerun/kegg.Rds")
  go %>% write_rds("data/input_processed/humann_rerun/go.Rds")
  l4ec %>% write_rds("data/input_processed/humann_rerun/l4ec.Rds")
  
  n_humann_annotations <-
    lapply(c("genefamilies", "pathabundance", "pathcoverage", "eggnog", "go", "ko", "level4ec", "pfam"), function(dat) {
      read_tsv(paste0("data/input_preprocessed/metagenome/humann_rerun/",dat,"_stats.tsv"), col_types = cols()) %>% 
        rename(sample_id = 1, n_annotations = 2) %>% 
        mutate(sample_id = str_extract(sample_id, "^S_[:digit:]*")) %>% 
        mutate(annotation_type = dat)
    }) %>% 
    bind_rows() %>% 
    pivot_wider(names_from = annotation_type, values_from = n_annotations)
  
  n_humann_annotations %>% 
    write_tsv("data/input_processed/humann_rerun/n_humann_annotations.tsv")
  
  # tmp <-
  #   pathways %>% 
  #   pivot_longer(-sample_id, values_to = "orig_pw") %>% 
  #   left_join(pathways_2 %>% 
  #               pivot_longer(-sample_id, values_to = "new_pw"))
  # 
  # tmp %>% 
  #   filter(name %in% tmp$name[!is.na(tmp$new_pw)][1:16]) %>% 
  #   ggplot(aes(x = orig_pw, y = new_pw)) +
  #   geom_point() +
  #   facet_wrap(~name) +
  #   geom_abline(slope = 1, intercept = 0)
  
  
  # sample_meta -------------------------------------------------------------
  
  sample_meta %>% 
    mutate(exclude_seq_coverage = Total_Bases_QC_ATLAS >= 1e9) %>% 
    filter(Total_Bases_QC_ATLAS >= 1e9,
           Prøvetype %in% "Baseline") %>% 
    write_rds("data/input_processed/sample_data.Rds")
  
  sample_meta %>% 
    mutate(exclude_seq_coverage = Total_Bases_QC_ATLAS >= 1e9) %>% 
    write_rds("data/input_processed/sample_meta.Rds")
  
  
  # screening_data ----------------------------------------------------------
  
  screening_data %>% 
    select(deltaker_id) %>% 
    left_join(lesions_data %>% select(deltaker_id, pathology_processed, side_processed), by = "deltaker_id") %>% 
    group_by(deltaker_id) %>% 
    summarize(proximal_non_advanced_adenoma = ifelse(any(str_detect(pathology_processed, "^5a") & side_processed %in% "Proximal"), 1, 0),
              distal_non_advanced_adenoma = ifelse(any(str_detect(pathology_processed, "^5a") & side_processed %in% "Distal"), 1, 0),
              .groups = "drop") %>% 
    right_join(screening_data, by = "deltaker_id") %>% 
    mutate(proximal_acn = case_when(proximal_crc %in% 1 ~ 1,
                                    proximal_advanced_adenoma %in% 1 ~ 1,
                                    proximal_advanced_serrated %in% 1 ~ 1,
                                    TRUE ~ 0),
           distal_acn = case_when(distal_crc %in% 1 ~ 1,
                                  distal_advanced_adenoma %in% 1 ~ 1,
                                  distal_advanced_serrated %in% 1 ~ 1,
                                  TRUE ~ 0)) %>%
    left_join(lesions_size_by_participant, by = join_by(deltaker_id)) %>% 
    left_join(lesions_size_data %>% 
                left_join(lesions_data %>% 
                            select(deltaker_id, lesionno, side_processed), 
                          by = join_by(deltaker_id, lesionno)) %>% 
                group_by(deltaker_id) %>% 
                summarize(sum_lesion_diameter = sum(diameter),
                          sum_lesion_diameter_proximal = ifelse(any(side_processed == "Proximal"), sum(diameter[ side_processed == "Proximal"]), 0),
                          sum_lesion_diameter_distal = ifelse(any(side_processed == "Distal"), sum(diameter[ side_processed == "Distal"]), 0),
                          .groups = "drop"), by = join_by(deltaker_id)) %>% 
    mutate(across(starts_with("sum_lesion_diameter"), function(x) {
      x[ is.na(x)] <- 0
      x
    })) %>% 
    mutate(final_result_cat_neg = case_when(detect_worthy_lesions == "negative" ~ "Negative",
                                            final_result == "5b. >= 3 Non-advanced adenomas" ~ ">= 3 Non-advanced adenomas",
                                            final_result == "5d. Advanced serrated" ~ "Advanced serrated",
                                            final_result == "5c. Advanced adenoma" ~ "Advanced adenoma",
                                            final_result == "6. Cancer" ~ "Cancer") %>% 
             factor(levels = c("Negative", ">= 3 Non-advanced adenomas", "Advanced serrated", "Advanced adenoma", "Cancer"))) %>%
    mutate(across(starts_with("colo_") & !starts_with("colo_count"), 
                  .fns = function(x) x %>% factor(labels = c("No", "Yes")))) %>% 
    write_rds("data/input_processed/screening_data.Rds")
  
  prior_samples %>% 
    mutate(first_round = case_when(prior_samples %in% 0 ~ "Yes",
                                   TRUE ~ "No")) %>% 
    write_tsv("data/input_processed/prior_samples.tsv")
  
  
  screening_data %>% 
    select(deltaker_id, final_result) %>% 
    left_join(lesions_data, by = "deltaker_id") %>% 
    ## Definition from Randel et al.
    mutate(adv_serrated_lesion = case_when(str_detect(lesiontype, "(serrated|Hyperplastisk)") & 
                                             (diameter >= 10 | 
                                                str_detect(dysplasiadiff, "^(Lett|Grov).*dysplasi$")) ~ TRUE,
                                           TRUE ~ FALSE)) %>% 
    left_join(lesions_size_data %>% select(-diameter), by = join_by(deltaker_id, lesionno)) %>% 
    mutate(lesion_diag = case_when(!adv_serrated_lesion ~ pathology_processed,
                                   adv_serrated_lesion ~ "5d. Advanced serrated lesion")) %>% 
    select(deltaker_id, locX, locY, lesion_diag, segment, segment_id, side_processed, diameter, polyp_volume = Polyp_vol) %>% 
    write_tsv("data/input_processed/lesion_locs.tsv")
  
    # count(final_result, fin_res) %>% 
    # pivot_wider(names_from = fin_res, values_from = n, values_fill = 0)
    # group_by(deltaker_id) %>% 
    # mutate(any_adv_serr = any(adv_serrated_lesion)) %>% 
    # ungroup() %>% 
    # filter(str_detect(final_result, "^5d."), !any_adv_serr) %>% 
    # View()
    # count(deltaker_id, adv_serrated_lesion, final_result) %>%
    # pivot_wider(names_from = adv_serrated_lesion, values_from = n, values_fill = 0) %>%
    # View()
    
  
  # lifestyle ---------------------------------------------------------------
  
  lifestyle_data %>%
    write_rds("data/input_processed/lifestyle_data.Rds")
  
  # diet --------------------------------------------------------------------
  
  diet_data %>% 
    left_join(screening_data %>% select(id, kjonn)) %>% 
    filter(exclude_reserved %in% 0) %>% 
    filter(exclude_nocolo %in% 0) %>% 
    skim(kjonn, exclude_nocolo, exclude_reserved, exclude_energyhigh, exclude_energylow, exclude_FFQquality_low)
  
  ## Summarize exclusion
  diet_data %>% 
    filter(exclude_nocolo %in% 0,
           exclude_reserved %in% 0) %>% ## FFQ for those meeting general inclusion criteria = 1562 participants
    summarize(FFQ_qual_low = sum(exclude_FFQquality_low %in% 1),
              energy_low = sum(exclude_energylow %in% 1),
              energy_high = sum(exclude_energyhigh %in% 1),
              remaining = sum(exclude_FFQquality_low %in% 0 & exclude_energylow %in% 0 & exclude_energyhigh %in% 0)) ## Number excluded by reason for participants with FFQs
  
  ## Exclude those not meeting criteria:
  diet_data_proc <- diet_data %>%
    filter(exclude_nocolo %in% 0,
           exclude_reserved %in% 0,
           exclude_FFQquality_low %in% 0,
           exclude_energyhigh %in% 0,
           exclude_energylow %in% 0) %>% 
    mutate(BMI_cat = case_when(!is.na(BMI_cat3) ~ as.character(BMI_cat3),
                               is.na(BMI_cat3) ~ "Missing")) %>% 
    mutate(BMI_cat = factor(BMI_cat, levels = c("Normal weight", "Overweight", "Obese", "Missing")))
  
  diet_data_proc %>% 
    write_rds("data/input_processed/diet_data.Rds")
  
  
  # diet wcrf ---------------------------------------------------------------
  
  diet_wcrf %>% 
    write_rds("data/input_processed/diet_wcrf.Rds")
  
  
  # diet energy adjusted ----------------------------------------------------
  
  diet_data_proc %>% 
    select(id) %>% 
    left_join(diet_ea) %>% 
    write_rds("data/input_processed/diet_ea.Rds")
  

  # antibiotics -------------------------------------------------------------

  antibiotics_filtered_from_Sandra <- 
    read_rds("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/papers/Microbiome_stability/datasets/metadata/antibiotics/Antibiotics_filtered_fromSandra.rds") %>% 
    tibble() %>% 
    select(-any_of(names(screening_data %>% select(-deltaker_id))))
  
  get_antibiotics <- function() {
    tmp <-
      antibiotics_filtered_from_Sandra %>% 
      left_join(screening_data %>% select(deltaker_id, fobt_investigation_date), by = "deltaker_id") %>% 
      mutate(months_since_sampling_1 = round((utleveringsdato - fobt_investigation_date)/(365.25/12), 1)) %>% 
      group_by(deltaker_id) %>% 
      summarize(first = min(months_since_sampling_1), 
                last = max(months_since_sampling_1),
                tot_prescriptions = n(),
                tot_DDD = sum(Utlevering_DDD),
                n_diff = length(unique(Legemiddel_ATCkode_Niva5)),
                DDD_year_before_fobt = sum(Utlevering_DDD[ months_since_sampling_1>-12 & months_since_sampling_1<0]),
                prescribed_4mo_before_fobt = last > -4, 
                .groups = "drop")
      
    tmp %>% 
      bind_rows(screening_data %>% 
                  select(deltaker_id) %>% 
                  filter(!deltaker_id %in% tmp$deltaker_id) %>% 
                  mutate(tot_prescriptions = 0,
                         tot_DDD = 0,
                         n_diff = 0,
                         DDD_year_before_fobt = 0,
                         prescribed_4mo_before_fobt = FALSE)) %>% 
      left_join(screening_data %>% 
                select(deltaker_id, id) %>% 
                left_join(lifestyle_data %>% 
                            select(id, Antibiotics), by = "id") %>% 
                select(deltaker_id, self_rep_antibiotics = Antibiotics), by = "deltaker_id") %>% 
      mutate(antibiotics_reg_quest_comb = case_when(prescribed_4mo_before_fobt | self_rep_antibiotics %in% "Yes" ~ "Yes",
                                                     TRUE ~ "No"))
  }
  
  antibiotics_summary <- get_antibiotics()
  
  antibiotics_summary %>% 
    write_tsv("data/input_processed/antibiotics_summary_data.tsv")
  

  # PPI ---------------------------------------------------------------------

  lmr_crcbiome <- read_csv2("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/62_161_uttrekksfil_2023-10-26T16.24.21.csv")
  lmr_key_crcbiome <- read_csv("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/crcbiome_lmr_nokkel_2022-06-15.csv")
  
  # lmr_crcbiome %>% 
  #   count(Legemiddel_ATCkode_Niva5) %>% 
  #   filter(str_detect(Legemiddel_ATCkode_Niva5, "^A0(1|2)"))  
  
  get_PPI <- function(only_PPI = FALSE) {
    
    if (!only_PPI) {
      re <- "^A02(A|B)"
    } else {
      re <- "^A02BC"
    }
    
    tmp <- 
      lmr_crcbiome %>% 
      ## Will include both Antacids and PPI prescriptions
      filter(str_detect(Legemiddel_ATCkode_Niva5, re)) %>% 
      mutate(Utlevering_ArManed = as_date(paste0(Utlevering_ArManed, "-15"), format = "%Y-%m-%d")) %>% 
      left_join(lmr_key_crcbiome, by = "lmr_kobling_id") %>% 
      left_join(screening_data %>% select(deltaker_id, fobt_investigation_date), by = "deltaker_id") %>% 
      filter(fobt_investigation_date > Utlevering_ArManed) %>% 
      group_by(deltaker_id) %>% 
      summarize(days_since_last = min(difftime(fobt_investigation_date, Utlevering_ArManed, units = "days")) %>% as.integer(),
                n_prescriptions_1y = sum(difftime(fobt_investigation_date, Utlevering_ArManed, units = "days") < 365.25),
                n_prescriptions_5y = sum(difftime(fobt_investigation_date, Utlevering_ArManed, units = "days") < 5*365.25),
                ddd_1y = sum(Utlevering_DDD[ difftime(fobt_investigation_date, Utlevering_ArManed, units = "days") < 365.25]),
                ddd_5y = sum(Utlevering_DDD[ difftime(fobt_investigation_date, Utlevering_ArManed, units = "days") < 5*365.25]),
                .groups = "drop") %>% 
      mutate(prescribed_4mo_before_fobt = factor(as.integer(days_since_last < 30.5*4), labels = c("No", "Yes")),
             any_PPI = "Yes") 
    
    tmp %>% 
      bind_rows(screening_data %>% 
                  select(deltaker_id) %>% 
                  filter(!deltaker_id %in% tmp$deltaker_id) %>% 
                  mutate(n_prescriptions_1y = 0,
                         n_prescriptions_5y = 0,
                         ddd_1y = 0,
                         ddd_5y = 0,
                         prescribed_4mo_before_fobt = "No",
                         any_PPI = "No")) %>% 
      left_join(screening_data %>% 
                  select(deltaker_id, id) %>% 
                  left_join(lifestyle_data %>% 
                              select(id, Antacids), by = "id") %>% 
                  select(deltaker_id, self_rep_antacids = Antacids), by = "deltaker_id") %>% 
      mutate(PPI_antacids_reg_quest_comb = case_when(prescribed_4mo_before_fobt %in% "Yes" | self_rep_antacids %in% "Yes" ~ "Yes",
                                                     TRUE ~ "No"))  
  }
  
  PPI_dat <- get_PPI()
  get_PPI(only_PPI = TRUE) %>% write_tsv("data/misc/PPI_use.tsv")
  
  rm(get_PPI)
    
  PPI_dat %>% 
    write_tsv("data/input_processed/PPI_data.tsv")
  
  
  lmr_key_crcbiome <- read_csv("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/crcbiome_lmr_nokkel_2022-06-15.csv")
  
  lmr_crcbiome <- 
    read_csv2("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/62_161_uttrekksfil_2023-10-26T16.24.21.csv") %>% 
    mutate(Utlevering_ArManed = as.Date(paste0(Utlevering_ArManed, "-15"))) %>% 
    left_join(lmr_key_crcbiome, by = "lmr_kobling_id") %>% 
    left_join(screening_data %>% select(deltaker_id, fobt_investigation_date, mb_invitert_dato), by = "deltaker_id") %>% 
    select(id = deltaker_id, everything()) %>% 
    mutate(months_since_invitation = as.double(round((Utlevering_ArManed - fobt_investigation_date)/(365.25/12), 1)))
  
  lmr_crcbiome_bcsn_invit <-
    read_csv2("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/62_161_uttrekksfil_2023-10-26T16.24.21.csv", 
              col_types = cols()) %>% 
    mutate(Utlevering_ArManed = as.Date(paste0(Utlevering_ArManed, "-15"))) %>% 
    left_join(lmr_key_crcbiome, by = "lmr_kobling_id") %>% 
    left_join(first_round_invitation %>% 
                select(deltaker_id, postlagt_dato, age_invitation),
              by = join_by(deltaker_id)) %>% 
    select(id = deltaker_id, everything()) %>% 
    mutate(months_since_invitation = as.double(round((Utlevering_ArManed - postlagt_dato)/(365.25/12), 1)))
  
  ## Mock invitation dates for population control group
  pop_mock_invitations_crcbiome_matched <- read_tsv("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/mock_invitations_pop_control_crcbiome_matched.tsv")
  pop_mock_invitations_screen_rep <- read_tsv("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/mock_invitations_pop_control_screen_rep.tsv")
  
  
  lmr_pop_screen_rep <- 
    read_csv2("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/62_169_uttrekksfil_2023-11-02T11.22.29.csv") %>% 
    mutate(Utlevering_ArManed = as.Date(paste0(Utlevering_ArManed, "-15"))) %>% 
    left_join(pop_mock_invitations_screen_rep %>% select(-c(Pasient_Fodselsar, Pasient_Kjonn_Verdi)), by = "PID") %>% 
    select(id = PID, everything()) %>% 
    mutate(months_since_invitation_scr_rep = as.integer(round((Utlevering_ArManed - mock_invitation_date_screening_rep)/(365.25/12), 1))) #%>% 
  
  lmr_pop_crcbiome_match <- 
    read_csv2("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/62_169_uttrekksfil_2023-11-02T11.22.29.csv") %>% 
    mutate(Utlevering_ArManed = as.Date(paste0(Utlevering_ArManed, "-15"))) %>% 
    inner_join(pop_mock_invitations_crcbiome_matched %>% select(-c(Pasient_Fodselsar, Pasient_Kjonn_Verdi, crcbiome_match)), by = "PID", relationship = "many-to-many") %>% 
    select(id = unique_resampling_id, -PID, everything()) %>% 
    mutate(months_since_invitation_crcbiome_match = as.integer(round((Utlevering_ArManed - mock_invitation_date_crcbiome_match)/(365.25/12), 1))) #%>% 
  
  ## Invitation dates and screening center for BCSN participants
  bcsn_invitation_data <- 
    read_csv2("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/crcbiome_lmr_kontrollgruppe_data_tillegg_avid.csv") %>% 
    left_join(read_csv("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/crcbiome_lmr_kontrollgruppe_fylke_avid.csv"), by = "PID")
  
  ## LMR data for BCSN participants
  lmr_BSCN <- 
    read_csv2("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/62_162_uttrekksfil_2023-10-26T16.25.53.csv") %>% 
    mutate(Utlevering_ArManed = as.Date(paste0(Utlevering_ArManed, "-15"))) %>%
    left_join(bcsn_invitation_data, by = "PID") %>%
    select(id = PID, everything()) %>%
    mutate(months_since_invitation = as.integer(round((Utlevering_ArManed - postlagt_dato1)/(365.25/12), 1)))
  
  
  summarize_lmr_comorbid_data <- function(dataset, n_months, id_dataset) {
    dataset %>% 
      group_by(id) %>% 
      filter(months_since_invitation <= 0,
             months_since_invitation > -n_months) %>% 
      summarize(
        ## Diabetes
        # diab_1_pres = any(str_detect(Legemiddel_ATCkode_Niva5, "^A10")),
        diab_2_pres = sum(str_detect(Legemiddel_ATCkode_Niva5, "^A10")) > 1,
        ## CVD: B01 or any C
        cvd_pres = any(str_detect(Legemiddel_ATCkode_Niva5, "^(B01|C)")),
        ## COPD
        copd_pres = any(str_detect(Legemiddel_ATCkode_Niva5, "^(R03AC|R03AK|R03AL|R03BB|R03DA|R03DX|R03DC|R03BA)")),
        ## Antibiotics: J01
        antibiotics_pres = any(str_detect(Legemiddel_ATCkode_Niva5, "^(J01)")),
        ## Antibiotics 2: J01, A07AA P01AB
        antibiotics_pres_2 = any(str_detect(Legemiddel_ATCkode_Niva5, "^(J01|A07AA|P01AB)")),
        .groups = "drop") %>% 
      right_join(id_dataset, by = "id") %>% 
      mutate(across(.cols = -id, .fns = function(x) ifelse(is.na(x), FALSE, x)))
  }
  
  
  crcbiome_comorbid_prescriptions <- 
    lmr_crcbiome %>% 
    summarize_lmr_comorbid_data(n_months = 12, id_dataset = screening_data %>% select(id = deltaker_id)) %>% 
    pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
    mutate(within = "12 months") %>% 
    bind_rows(lmr_crcbiome %>% 
                summarize_lmr_comorbid_data(n_months = 24, id_dataset = screening_data %>% select(id = deltaker_id)) %>% 
                pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
                mutate(within = "24 months")) %>% 
    left_join(screening_data %>% 
                select(id = deltaker_id, sex = kjonn, fobt_investigation_date) %>% 
                left_join(lmr_crcbiome %>% 
                            select(id, Pasient_Fodselsar) %>% 
                            distinct(), by = "id") %>% 
                mutate(age = year(fobt_investigation_date)-Pasient_Fodselsar) %>% 
                select(id, age, sex), 
              by = "id")
  
  crcbiome_screening_rep_comorbid_prescriptions <-
    lmr_crcbiome_bcsn_invit  %>%
    # filter(crcbiome_consent %in% 1) %>% 
    summarize_lmr_comorbid_data(n_months = 12, id_dataset = screening_data %>% select(id = deltaker_id)) %>% 
    pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
    mutate(within = "12 months") %>% 
    bind_rows(lmr_crcbiome_bcsn_invit %>% 
                summarize_lmr_comorbid_data(n_months = 24, id_dataset = screening_data %>% select(id = deltaker_id)) %>% 
                pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
                mutate(within = "24 months")) %>% 
    left_join(first_round_invitation %>% 
                select(id = deltaker_id, sex = kjonn, age = age_invitation) %>% 
                mutate(sex = case_when(sex == "F" ~ "Female",
                                       sex == "M" ~ "Male")), 
              by = "id")
  
  bcsn_no_crcbiome_comorbid_prescriptions <-
    lmr_BSCN  %>%
    filter(crcbiome_consent %in% 0) %>%
    summarize_lmr_comorbid_data(n_months = 12, id_dataset = bcsn_invitation_data %>% filter(crcbiome_consent %in% 0) %>% select(id = PID)) %>% 
    pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
    mutate(within = "12 months") %>% 
    bind_rows(lmr_BSCN %>% 
                summarize_lmr_comorbid_data(n_months = 24, id_dataset = bcsn_invitation_data %>% filter(crcbiome_consent %in% 0) %>% select(id = PID)) %>% 
                pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
                mutate(within = "24 months")) %>% 
    left_join(bcsn_invitation_data %>% 
                select(id = PID, everything()) %>% 
                mutate(sex = case_when(kjonn == "F" ~ "Female",
                                       kjonn == "M" ~ "Male"),
                       arm = case_when(arm == "F" ~ "FIT",
                                       arm == "S" ~ "Sigmoidoscopy"),
                       age = year(postlagt_dato1) - fodselsaar) %>% 
                select(-c(kjonn, fodselsaar, crcbiome_consent, crcbiome_colonoscopy, postlagt_dato1)),
              # select(id, age), 
              by = join_by(id))
  
  bcsn_comorbid_prescriptions <-
    lmr_BSCN  %>%
    filter(crcbiome_consent %in% c(0,1)) %>%
    summarize_lmr_comorbid_data(n_months = 12, id_dataset = bcsn_invitation_data %>% filter(crcbiome_consent %in% c(0,1)) %>% select(id = PID)) %>% 
    pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
    mutate(within = "12 months") %>% 
    bind_rows(lmr_BSCN %>% 
                summarize_lmr_comorbid_data(n_months = 24, id_dataset = bcsn_invitation_data %>% filter(crcbiome_consent %in% c(0,1)) %>% select(id = PID)) %>% 
                pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
                mutate(within = "24 months")) %>% 
    left_join(bcsn_invitation_data %>% 
                select(id = PID, everything()) %>% 
                mutate(sex = case_when(kjonn == "F" ~ "Female",
                                       kjonn == "M" ~ "Male"),
                       arm = case_when(arm == "F" ~ "FIT",
                                       arm == "S" ~ "Sigmoidoscopy"),
                       age = year(postlagt_dato1) - fodselsaar) %>% 
                select(-c(kjonn, fodselsaar, crcbiome_consent, crcbiome_colonoscopy, postlagt_dato1)),
              # select(id, age), 
              by = join_by(id))

  # bcsn_comorbid_prescriptions <-
  #   lmr_BSCN  %>%
  #   summarize_lmr_comorbid_data(n_months = 12, id_dataset = bcsn_invitation_data %>% select(id = PID)) %>% 
  #   pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
  #   mutate(within = "12 months") %>% 
  #   bind_rows(lmr_BSCN %>% 
  #               summarize_lmr_comorbid_data(n_months = 24, id_dataset = bcsn_invitation_data %>% select(id = PID)) %>% 
  #               pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
  #               mutate(within = "24 months")) %>% 
  #   left_join(bcsn_invitation_data %>% 
  #               select(id = PID, everything()) %>% 
  #               left_join(lmr_BSCN %>% 
  #                           select(id, Pasient_Fodselsar, sex = Pasient_Kjonn_Verdi) %>% 
  #                           distinct(), by = "id") %>% 
  #               mutate(age = year(postlagt_dato1)-Pasient_Fodselsar) %>% 
  #               select(id, age, sex),
  #             by = "id") %>% 
  #   mutate(sex = case_when(sex == 1 ~ "Male",
  #                          sex == 2 ~ "Female")) %>% 
  #   select(-age) %>% 
  #   left_join(read_csv("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/lmr/crcbiome_lmr_kontrollgruppe_fylke_avid.csv") %>%
  #               left_join(bcsn_invitation_data, by = "PID") %>% 
  #               select(id = PID, everything()) %>% 
  #               mutate(sex = case_when(kjonn == "F" ~ "Female",
  #                                      kjonn == "M" ~ "Male"),
  #                      arm = case_when(arm == "F" ~ "FIT",
  #                                      arm == "S" ~ "Sigmoidoscopy"),
  #                      age = year(postlagt_dato1) - fodselsaar) %>% 
  #               select(-c(kjonn, fodselsaar, crcbiome_consent, crcbiome_colonoscopy, postlagt_dato1)),
  #               # select(id, age), 
  #             by = join_by(id)) %>% 
  #     select(-c(fodselsaar, crcbiome_consent, crcbiome_colonoscopy))
  
  ## Use date for screening representative population
  pop_comorbid_prescriptions_scr_rep <-
    lmr_pop_screen_rep %>% 
    rename(months_since_invitation = months_since_invitation_scr_rep) %>% 
    summarize_lmr_comorbid_data(n_months = 12, id_dataset = pop_mock_invitations_screen_rep %>% select(id = PID)) %>% 
    pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
    mutate(within = "12 months") %>% 
    bind_rows(lmr_pop_screen_rep %>% 
                rename(months_since_invitation = months_since_invitation_scr_rep) %>% 
                summarize_lmr_comorbid_data(n_months = 24, id_dataset = pop_mock_invitations_screen_rep %>% select(id = PID)) %>% 
                pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
                mutate(within = "24 months")) %>% 
    left_join(pop_mock_invitations_screen_rep %>% 
                select(id = PID, sex = Pasient_Kjonn_Verdi, everything()) %>% 
                mutate(age = year(mock_invitation_date_screening_rep)-Pasient_Fodselsar) %>% 
                select(id, age, sex),
              by = "id") %>% 
    mutate(sex = case_when(sex == 1 ~ "Male",
                           sex == 2 ~ "Female"))
  
  ## Use date for crcbiome representative population
  # pop_comorbid_prescriptions_crcbiome_crcbiome_match <-
  #   lmr_pop_crcbiome_match %>% 
  #   rename(months_since_invitation = months_since_invitation_crcbiome_match) %>% 
  #   summarize_lmr_comorbid_data(n_months = 12, id_dataset = pop_mock_invitations_crcbiome_matched %>% select(id = unique_resampling_id)) %>% 
  #   pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
  #   mutate(within = "12 months") %>% 
  #   bind_rows(lmr_pop_crcbiome_match %>% 
  #               rename(months_since_invitation = months_since_invitation_crcbiome_match) %>% 
  #               summarize_lmr_comorbid_data(n_months = 24, id_dataset = pop_mock_invitations %>% select(id = unique_resampling_id)) %>% 
  #               pivot_longer(-id, names_to = "drug_cat", values_to = "prescribed") %>% 
  #               mutate(within = "24 months")) %>% 
  #   left_join(pop_mock_invitations_crcbiome_matched %>% 
  #               select(id = unique_resampling_id, everything()) %>% 
  #               mutate(age = year(mock_invitation_date_crcbiome_match)-Pasient_Fodselsar) %>% 
  #               select(id, age, sex = Pasient_Kjonn_Verdi),
  #             by = "id") %>% 
  #   mutate(sex = case_when(sex == 1 ~ "Male",
  #                          sex == 2 ~ "Female"))
  
  
  crcbiome_screening_rep_comorbid_prescriptions %>% 
    write_tsv("data/lmr/crcbiome_screening_rep_comorbid_prescriptions.tsv")
  
  bcsn_no_crcbiome_comorbid_prescriptions %>% 
    write_tsv("data/lmr/bcsn_no_crcbiome_comorbid_prescriptions.tsv")
  
  crcbiome_comorbid_prescriptions %>% 
    write_tsv("data/lmr/crcbiome_comorbid_prescriptions.tsv")
  
  bcsn_comorbid_prescriptions %>% 
    write_tsv("data/lmr/bcsn_comorbid_prescriptions.tsv")
  
  pop_comorbid_prescriptions_scr_rep %>% 
    write_tsv("data/lmr/pop_comorbid_prescriptions_scr_rep.tsv")
  
  # pop_comorbid_prescriptions_crcbiome_crcbiome_match %>% 
  #   write_tsv("data/lmr/pop_comorbid_prescriptions_crcbiome_match.tsv")
  
  
  # metaphlan_species -------------------------------------------------------
  
  sample_meta %>% 
    filter(Total_Bases_QC_ATLAS >= 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id) %>% 
    left_join(abundance_data) %>% 
    ## Restrict dataset to those present in any samples
    select(1, which(apply(.[,-1]>0, 2, any))+1) %>% 
    write_rds("data/input_processed/metaphlan_species.Rds")
  
  sample_meta %>% 
    filter(Total_Bases_QC_ATLAS >= 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id) %>% 
    left_join(mp4_sgbs) %>% 
    select(1, which(apply(.[,-1]>0, 2, any))+1) %>% 
    write_rds("data/input_processed/mp4_sgbs.Rds")
  
  sample_meta %>% 
    filter(Total_Bases_QC_ATLAS >= 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id) %>% 
    left_join(mp4_species) %>% 
    select(1, which(apply(.[,-1]>0, 2, any))+1) %>% 
    write_rds("data/input_processed/mp4_species.Rds")
  
  sample_meta %>% 
    filter(Total_Bases_QC_ATLAS >= 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id) %>% 
    left_join(mp4_genera) %>% 
    select(1, which(apply(.[,-1]>0, 2, any))+1) %>% 
    write_rds("data/input_processed/mp4_genera.Rds")
  
  mp4_sgb_names %>% 
    # mutate(phylum = str_replace(phylum, "Firmicutes", "Bacillota")) %>% 
    # mutate(clade_name = str_replace(clade_name, "p__Firmicutes", "p__Bacillota")) %>% 
    write_tsv("data/input_processed/mp4_sgb_names.tsv")
  mp4_species_names %>% 
    # mutate(phylum = str_replace(phylum, "Firmicutes", "Bacillota")) %>% 
    # mutate(clade_name = str_replace(clade_name, "p__Firmicutes", "p__Bacillota")) %>% 
    write_tsv("data/input_processed/mp4_species_names.tsv")
  mp4_genera_names %>% 
    # mutate(phylum = str_replace(phylum, "Firmicutes", "Bacillota")) %>% 
    # mutate(clade_name = str_replace(clade_name, "p__Firmicutes", "p__Bacillota")) %>% 
    write_tsv("data/input_processed/mp4_genera_names.tsv")
  
  

  # colibactin --------------------------------------------------------------
  
  pks_mapping <- 
    read_tsv("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/masterprojects/Lisi_masterproject/pipeline3_read_mapping/results/coverage_pks.tsv",
             col_names = c("covered_bases", "coverage", "mean_coverage", "sample_id")) %>% 
    mutate(sample_id = str_replace(sample_id, "-", "_"))

  ecoli_mapping <- 
    read_tsv("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/masterprojects/Lisi_masterproject/pipeline3_read_mapping/results/coverage_ecoli.tsv",
             col_names = c("covered_bases", "coverage", "mean_coverage", "sample_id")) %>% 
    mutate(sample_id = str_replace(sample_id, "-", "_"))
  
  pks_mapping %>% 
    select(pks_covered_bases = covered_bases, sample_id) %>% 
    left_join(ecoli_mapping %>% 
                select(ecoli_covered_bases = covered_bases, sample_id), by = join_by(sample_id)) %>% 
    inner_join(sample_data %>% select(sample_id), by = join_by(sample_id)) %>% 
    mutate(pks_pos = pks_covered_bases > 0.1*max(pks_covered_bases)) %>% 
    mutate(pks_frac = pks_covered_bases/max(pks_covered_bases)) %>% 
    mutate(ecoli_frac = ecoli_covered_bases/max(ecoli_covered_bases)) %>% 
    mutate(pks_dominating = ecoli_frac < pks_frac + 0.15 &
             ecoli_frac > pks_frac-0.15) %>% 
    ggplot(aes(x = pks_covered_bases, y = ecoli_covered_bases, color = pks_dominating)) +
    facet_wrap(~pks_pos) +
    geom_point()
   
  pks_mapping %>% 
    inner_join(sample_data %>% select(sample_id), by = "sample_id") %>% 
    mutate(partial_pks = case_when(coverage < 97.5 & mean_coverage > 7 ~ "partial_pks",
                                   coverage < 10 & mean_coverage > 0.2  ~ "small_partial_pks",
                                   coverage < 5 & mean_coverage > 0.07  ~ "small_partial_pks",
                                   coverage > 97.5 ~ "complete",
                                   coverage < 10 ~ "not_detected",
                                   TRUE ~ "looks_normal")) %>% 
    select(sample_id, pks_dp = mean_coverage, partial_pks, covered_bases) %>%
    ## Use metaphlan data to determine if e coli is present
    inner_join(sample_data %>% 
                select(sample_id) %>% 
                left_join(mp4_sgbs, by = "sample_id") %>% 
                select(sample_id, any_of(mp4_sgb_names %>% 
                                           filter(str_detect(species, "Escherichia_coli")) %>% 
                                           pull(sgb))) %>% 
                rename("mp4_ecoli" = 2), by = "sample_id") %>% 
    left_join(ecoli_mapping %>% 
                select(ecoli_dp = mean_coverage, 
                       ecoli_covered_bases = covered_bases,
                       sample_id), by = join_by(sample_id)) %>% 
    mutate(ecoli_pos = ifelse(mp4_ecoli > 0, "present", "absent")) %>%
    mutate(pks_classification = case_when(ecoli_pos %in% "absent" ~ "genotoxic pks negative", 
                                          partial_pks %in% c("partial_pks", "small_partial_pks", "not_detected") ~ "genotoxic pks negative",
                                          partial_pks %in% c("complete", "looks_normal") ~ "genotoxic pks positive")) %>% 
    # select(sample_id, pks_classification) %>% 
    write_tsv("data/input_processed/pks_detection.tsv")
    
  
  # MAGs --------------------------------------------------------------------
  
  MAGs <- 
    sample_meta %>%
    filter(Total_Bases_QC_ATLAS > 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id, Total_Reads_QC_ATLAS) %>% 
    left_join(median_coverage_genomes %>% 
                mutate(sample_id = gsub("-", "_", sample_id))) %>% 
    pivot_longer(-c(sample_id, Total_Reads_QC_ATLAS)) %>% 
    group_by(name) %>% 
    mutate(any_obs = any(value > 0)) %>% 
    filter(any_obs) %>% 
    select(-any_obs) %>%
    ungroup() %>% 
    mutate(value = value/Total_Reads_QC_ATLAS*1e6) %>% 
    select(-Total_Reads_QC_ATLAS) %>% 
    pivot_wider(names_from = "name", values_from = value)
  
  ## Remove unwanted samples and adjust the read counts of MAGs whose median depth in a sample is 0 to 0.
  ## Also remove low quality mags
  raw_counts_mags_filtered <-
    raw_counts_genomes %>% 
    pivot_longer(-Sample, names_to = "sample_id", values_to = "readcounts") %>% 
    mutate(sample_id =gsub("-", "_", sample_id)) %>% 
    filter(sample_id %in% MAGs$sample_id) %>% 
    filter(Sample %in% names(MAGs)) %>% 
    left_join(MAGs %>% 
                pivot_longer(-sample_id, names_to = "Sample", values_to = "median_depth")) %>% 
    mutate(filtered_read_counts = case_when(median_depth == 0 ~ 0,
                                            TRUE ~ readcounts)) %>% 
    group_by(Sample) %>% 
    mutate(any_obs = any(filtered_read_counts > 0)) %>% 
    ungroup() %>% 
    filter(any_obs) %>% 
    select(-c(readcounts, median_depth, any_obs)) %>% 
    pivot_wider(names_from = Sample, values_from = filtered_read_counts)
  
  raw_counts_mags_filtered %>% 
    write_rds("data/input_processed/MAG_readcounts_filtered.Rds")
  
  ## Remove unwanted samples and no adjustment to 0.
  ## Also remove low quality mags
  raw_counts_mags <-
    raw_counts_genomes %>% 
    pivot_longer(-Sample, names_to = "sample_id", values_to = "readcounts") %>% 
    mutate(sample_id =gsub("-", "_", sample_id)) %>% 
    filter(sample_id %in% MAGs$sample_id) %>% 
    filter(Sample %in% names(MAGs)) %>% 
    group_by(Sample) %>% 
    mutate(any_obs = any(readcounts > 0)) %>% 
    ungroup() %>% 
    filter(any_obs) %>% 
    select(-c(any_obs)) %>% 
    pivot_wider(names_from = Sample, values_from = readcounts)
  
  raw_counts_mags %>% 
    write_rds("data/input_processed/MAG_readcounts.Rds")
  
  mags_taxa_id %>% 
    filter(MAG_id %in% (MAGs %>% names())) %>% 
    select(MAG_id, domain, phylum, clade, order, family, genus, species, taxonomy_id) %>% 
    write_tsv("data/input_processed/MAG_taxonomy.tsv")
  
  MAGs %>% 
    write_rds("data/input_processed/MAG_abundance.Rds")
  
  
  # MAG eggnog --------------------------------------------------------------
  
  MAG_eggnog_norm <- sample_meta %>%
    filter(grepl("S_", sample_id),
           Total_Bases_QC_ATLAS > 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id, Total_Reads_QC_ATLAS) %>% 
    left_join(MAG_eggnog) %>% 
    select(-"NA") %>% 
    pivot_longer(-c(sample_id, Total_Reads_QC_ATLAS)) %>% 
    group_by(name) %>% 
    mutate(any_obs = any(value > 0)) %>% 
    filter(any_obs) %>% 
    select(-any_obs) %>%
    ungroup() %>% 
    mutate(value = value/Total_Reads_QC_ATLAS*1e6) %>% 
    select(-Total_Reads_QC_ATLAS) %>% 
    pivot_wider(names_from = "name", values_from = value)
  
  MAG_eggnog_norm %>% 
    write_rds("data/input_processed/MAG_eggnog_abundance.Rds")
  
  # MAG COGs ----------------------------------------------------------------
  
  MAG_COG_norm <- sample_meta %>%
    filter(grepl("S_", sample_id),
           Total_Bases_QC_ATLAS > 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id, Total_Reads_QC_ATLAS) %>% 
    left_join(MAG_COG) %>% 
    pivot_longer(-c(sample_id, Total_Reads_QC_ATLAS)) %>% 
    group_by(name) %>% 
    mutate(any_obs = any(value > 0)) %>% 
    filter(any_obs) %>% 
    select(-any_obs) %>%
    ungroup() %>% 
    mutate(value = value/Total_Reads_QC_ATLAS*1e6) %>% 
    select(-Total_Reads_QC_ATLAS) %>% 
    pivot_wider(names_from = "name", values_from = value)
  
  MAG_COG_norm %>% 
    write_rds("data/input_processed/MAG_COG_abundance.Rds")
  
  # MAG GOs -----------------------------------------------------------------
  
  MAG_GO_norm <- sample_meta %>%
    filter(grepl("S_", sample_id),
           Total_Bases_QC_ATLAS > 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id, Total_Reads_QC_ATLAS) %>% 
    left_join(MAG_GO) %>% 
    select(-starts_with("NA")) %>% 
    pivot_longer(-c(sample_id, Total_Reads_QC_ATLAS)) %>% 
    group_by(name) %>% 
    mutate(any_obs = any(value > 0)) %>% 
    filter(any_obs) %>% 
    select(-any_obs) %>%
    ungroup() %>% 
    mutate(value = value/Total_Reads_QC_ATLAS*1e6) %>% 
    select(-Total_Reads_QC_ATLAS) %>% 
    pivot_wider(names_from = "name", values_from = value)
  
  MAG_GO_norm %>% 
    write_rds("data/input_processed/MAG_GO_abundance.Rds")
  
  # MAG KOs -----------------------------------------------------------------
  
  MAG_KO_norm <- sample_meta %>%
    filter(grepl("S_", sample_id),
           Total_Bases_QC_ATLAS > 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id, Total_Reads_QC_ATLAS) %>% 
    left_join(MAG_KO) %>% 
    select(-starts_with("NA")) %>% 
    pivot_longer(-c(sample_id, Total_Reads_QC_ATLAS)) %>% 
    group_by(name) %>% 
    mutate(any_obs = any(value > 0)) %>% 
    filter(any_obs) %>% 
    select(-any_obs) %>%
    ungroup() %>% 
    mutate(value = value/Total_Reads_QC_ATLAS*1e6) %>% 
    select(-Total_Reads_QC_ATLAS) %>% 
    pivot_wider(names_from = "name", values_from = value)
  
  MAG_KO_norm %>% 
    write_rds("data/input_processed/MAG_KO_abundance.Rds")
  
  # viral abundance ---------------------------------------------------------
  
  ## Viruses
  virus_abundance <- 
    sample_meta %>%
    filter(Total_Bases_QC_ATLAS > 1e9,
           Prøvetype %in% "Baseline") %>% 
    select(sample_id, Total_Reads_QC_ATLAS) %>% 
    left_join(viral_abundance %>%
                pivot_longer(-ID, names_to = "sample_id") %>% 
                mutate(sample_id = gsub("-", "_", sample_id)) %>% 
                pivot_wider(names_from = ID, values_from = value)) %>% 
    pivot_longer(-c(sample_id, Total_Reads_QC_ATLAS)) %>% 
    group_by(name) %>% 
    mutate(any_obs = any(value > 0)) %>% 
    filter(any_obs) %>% 
    select(-any_obs) %>%
    ungroup() %>% 
    mutate(value = value/Total_Reads_QC_ATLAS*1e6) %>% 
    select(-Total_Reads_QC_ATLAS) %>% 
    pivot_wider(names_from = "name", values_from = value)
  
  virus_abundance %>% 
    write_rds("data/input_processed/viral_abundance.Rds")
  
  
  
  
  ## define variables
  
  variables_tmp <- bind_rows(
    c("Colonoscopy result" = "final_result",
      "CRC-related findings" = "detect_worthy_lesions",
      "Outcome" = "final_result_cat_neg",
      "Sex" = "kjonn", 
      "Age" = "age_invitation", 
      "Screening center" = "senter", 
      "FIT value" = "fobt_verdi", 
      "Programmatic screening round" = "runde",
      "Hemorrhoids" = "colo_hemorrhoids",
      "Diverticulitis" = "colo_diverticulitis",
      "IBD" = "colo_IBD",
      "BMI" = "BMI",
      "Height" = "Høyde", 
      "Physical activity" = "PhysAct_Score",
      "Smoking" = "Smoking", 
      "Snus" = "Snus",
      "Education" = "Utdanning", 
      "Marital status" = "Sivilstatus_cat2",
      "Employment status" = "Arbeid_lump",
      "National background" = "Nasj_cat2",
      "WCRF" = "wcrf_index_main",
      "Previous negative FIT" = "first_round",
      "Family history of CRC" = "Tarmkreft_Familie") %>% 
      enframe(name = "var_name", value = "var_id") %>% 
      mutate(dataset = case_when(var_name %in% c("BMI", "Høyde") ~ "Physical traits",
                                 var_name %in% c("Physical activity", "Smoking", "Snus", "WCRF") ~ "Lifestyle",
                                 var_name %in% c("final_result", "detect_worthy_lesions", "final_result_cat_neg", 
                                                 "fobt_verdi", "runde", "first_round",
                                                 "colo_hemorrhoids", "colo_diverticulitis", "colo_IBD") ~ "Screening-related",
                                 TRUE ~ "Demography")),
    
    c("Energy (kcal/day)" = "Energi_kcal", 
      "Red meat" = "ROKJOT_nonproc",
      "Vegetables" = "GRSAK_modified",
      "Fruits and berries" = "FRUKTB_modified",
      "Processed meat" = "PROKJOT_red",
      "Fish" = "FISK_modified",
      "Milk" = "MILK",
      "Cheeses" = "OST",
      "Sweet drinks" = "SABR_S",
      "Calcium" = "Ca",
      "Alcohol" = "Alko",
      "Fiber" = "Fiber_resid") %>% 
      enframe(name = "var_name", value = "var_id") %>% 
      mutate(dataset = "Diet"),
      
    c("Antibiotics use" = "antibiotics_reg_quest_comb",
      "Antacid/PPI use" = "PPI_antacids_reg_quest_comb",
      "Appendix removed" = "Blindtarm", 
      "Bowel disorder" = "Bowel_disorder_merged",
      "Self reported gastrointestinal disorders" = "bowel_disorder_merged_2") %>% 
      enframe(name = "var_name", value = "var_id") %>% 
      mutate(dataset = "Medical"),
    
    c("Abdominal pain" = "symtomduration_smerter_cat",
      "Rectal bleeding" = "symtomduration_synlig_blod_i_avforingen_cat",
      "Other symptoms" = "symtomduration_annet_se_fritekst_cat", 
      "Alternating bowel habits" = "symtomduration_vekslende_avforing_cat",
      "Change in bowel habits" = "symtomduration_endrede_avforingsvaner_cat",
      "Diarrhea" = "symtomduration_diare_cat",
      "Bloating" = "symtomduration_luftplager_cat",
      "Constipation" = "symtomduration_obstipasjon_cat",
      "Symptom burden" = "symptoms_overall_cat",
      "GI-related symptoms" = "symptoms_GI_related_cat") %>% 
      enframe(name = "var_name", value = "var_id") %>% 
      mutate(dataset = "Symptoms")) %>%
    
    mutate(lididem_crc = case_when(var_id %in% c("Energi_kcal", "Fiber_resid", "GRSAK_modified", 
                                                 "FRUKTB_modified", "ROKJOT_nonproc", "PROKJOT_red", 
                                                 "FISK_modified", "MILK", "OST", "SABR_S", 
                                                 "Ca", "Alko", "BMI", "Smoking", 
                                                 "PhysAct_Score", "kjonn", "age_invitation", "Høyde") ~ TRUE,
                                   TRUE ~ FALSE))
  
  meta_dat <- 
    read_rds("data/input_processed/sample_meta.Rds") %>% 
    filter(Prøvetype %in% "Baseline",
           Total_Bases_QC_ATLAS >= 1e9) %>% 
    select(deltaker_id) %>% 
    left_join(read_rds("data/input_processed/screening_data.Rds") %>% select(any_of(variables_tmp %>% pull(var_id)), deltaker_id, id), by = "deltaker_id") %>% 
    left_join(read_rds("data/input_processed/diet_data.Rds") %>% select(any_of(variables_tmp %>% pull(var_id)), id), by = "id") %>% 
    left_join(read_rds("data/input_processed/diet_ea.Rds") %>% select(-Energi_kcal) %>% select(any_of(variables_tmp %>% pull(var_id)), id), by = "id") %>% 
    left_join(read_rds("data/input_processed/diet_wcrf.Rds") %>% select(any_of(variables_tmp %>% pull(var_id)), id), by = "id") %>% 
    left_join(read_rds("data/input_processed/lifestyle_data.Rds") %>% select(any_of(variables_tmp %>% pull(var_id)), id), by = "id") %>% 
    left_join(read_tsv("data/input_processed/prior_samples.tsv", col_types = cols()) %>% select(any_of(variables_tmp %>% pull(var_id)), deltaker_id), by = "deltaker_id") %>% 
    left_join(read_tsv("data/input_processed/antibiotics_summary_data.tsv", col_types = cols()) %>% select(any_of(variables_tmp %>% pull(var_id)), deltaker_id), by = "deltaker_id") %>% 
    left_join(read_tsv("data/input_processed/PPI_data.tsv", col_types = cols()) %>% select(any_of(variables_tmp %>% pull(var_id)), deltaker_id), by = "deltaker_id") %>% 
    select(-id)
  
  variables <- 
    variables_tmp %>% 
    mutate(ref = case_when(var_id %in% "kjonn" ~ "Female",
                           var_id %in% "runde" ~ "2nd round",
                           var_id %in% "final_result" ~ "1. Negative",
                           var_id %in% "detect_worthy_lesions" ~ "negative",
                           var_id %in% "final_result_cat_neg" ~ "Negative",
                           var_id %in% "senter" ~ "Moss",
                           # var_id %in% "Antibiotics" ~ "No",
                           # var_id %in% "Antacids" ~ "No",
                           var_id %in% "Blindtarm" ~ "No",
                           var_id %in% "Smoking" ~ "Non smoker",
                           var_id %in% "Snus" ~ "Non snuser",
                           var_id %in% "Utdanning" ~ "Primary school",
                           var_id %in% "Sivilstatus_cat2" ~ "Married/cohabiting",
                           var_id %in% "Arbeid_lump" ~ "Employed",
                           var_id %in% "Nasj_cat2" ~ "Native",
                           var_id %in% "Bowel_disorder_merged" ~ "No bowel disease",
                           var_id %in% "bowel_disorder_merged_2" ~ "No",
                           var_id %in% "first_round" ~ "No",
                           var_id %in% "Tarmkreft_Familie" ~ "No",
                           var_id %in% "colo_diverticulitis" ~ "No",
                           var_id %in% "colo_hemorrhoids" ~ "No",
                           var_id %in% "colo_IBD" ~ "No",
                           str_detect(var_id, "symtom") ~ "No",
                           str_detect(var_id, "symptoms") ~ "No symptoms",
                           str_detect(var_id, "reg_quest") ~ "No")) %>% 
    mutate(var_id = factor(var_id, levels = names(meta_dat)[ names(meta_dat) %in% var_id]))
  
  variables %>% 
    write_rds("data/input_processed/metadata_variables.Rds")
  meta_dat %>% 
    write_rds("data/input_processed/metadata_selected_variables.Rds")
  meta_dat %>% 
    mutate(across(c(where(is.double), -deltaker_id), .fns = function(x) cat_func(x) %>% str_replace("negative", "low") %>% str_replace("positive", "high") %>% factor(levels = c("low", "mid", "high")))) %>% 
    write_rds("data/input_processed/metadata_selected_variables_cat.Rds")
  
  
  # variables <- 
  #   variables %>% 
  #   mutate(var_id = factor(var_id, levels = names(meta_dat)[ names(meta_dat) %in% var_id])) %>% 
  #   mutate(data_source = case_when(dataset %in% c("diet", "diet_ea", "diet_wcrf") ~ "FFQ",
  #                                  dataset %in% "lifestyle" ~ "LDQ",
  #                                  dataset %in% "screening" ~ "screening db"))
  
  lee_metadata <- 
    read_csv("data/input_preprocessed/Lee_2023/lee_metadata.csv") %>% 
    rename(sample_id = Run)
  
  # lee_metadata %>% 
  #   ggplot(aes(x = Bases, group = factor(AvgSpotLen), color = factor(AvgSpotLen))) +
  #   geom_density() +
  #   geom_vline(xintercept = 1e9)
  # 
  # lee_metadata %>% 
  #   mutate(under_1e9 = Bases < 1e9) %>% 
  #   count(under_1e9)
  
  lee_mp4 <- read_tsv("data/input_preprocessed/Lee_2023/profiles.tsv", skip = 1)
  
  lee_sgbs <-
    lee_mp4 %>% 
    filter(str_detect(clade_name, "t__")) %>% 
    mutate(clade_name = str_remove(clade_name, "^.*t__")) %>% 
    pivot_longer(-clade_name, names_to = "sample_id", values_to = "abundance") %>% 
    group_by(sample_id) %>% 
    mutate(abundance = abundance/sum(abundance)*100) %>% 
    ungroup() %>% 
    pivot_wider(names_from = clade_name, values_from = abundance)
  
  lee_species <-
    lee_mp4 %>% 
    filter(!str_detect(clade_name, "t__"),
           str_detect(clade_name, "s__")) %>% 
    mutate(clade_name = str_remove(clade_name, "\\|t__.*"),
           clade_name = str_remove(clade_name, ".*s__")) %>% 
    pivot_longer(-clade_name, names_to = "sample_id", values_to = "abundance") %>% 
    group_by(sample_id) %>% 
    mutate(abundance = abundance/sum(abundance)*100) %>% 
    ungroup() %>% 
    pivot_wider(names_from = clade_name, values_from = abundance)
  
  lee_sgb_names <- 
    lee_mp4 %>% 
    filter(str_detect(clade_name, "t__")) %>% 
    select(clade_name) %>% 
    separate(clade_name, into = c("kingdom", "phylum", "clade", "order", "family", "genus", "species", "sgb"), 
             sep = "\\|", 
             remove = FALSE) %>% 
    mutate(across(-clade_name, .fns = function(x) x %>% str_remove("[:alpha:]__")))
  
  lee_species_names <- 
    lee_mp4 %>% 
    filter(!str_detect(clade_name, "t__"),
           str_detect(clade_name, "s__")) %>% 
    select(clade_name) %>% 
    separate(clade_name, into = c("kingdom", "phylum", "clade", "order", "family", "genus", "species"), 
             sep = "\\|", 
             remove = FALSE) %>% 
    mutate(across(-clade_name, .fns = function(x) x %>% str_remove("[:alpha:]__")))
  
  lee_metadata %>% 
    filter(Bases >= 1e9) %>% 
    select(sample_id) %>% 
    left_join(lee_sgbs) %>% 
    select(1, which(apply(.[,-1]>0, 2, any))+1) %>% 
    write_rds("data/input_processed/lee_sgbs.Rds")
  
  lee_metadata %>% 
    filter(Bases >= 1e9) %>% 
    select(sample_id) %>% 
    left_join(lee_species) %>% 
    select(1, which(apply(.[,-1]>0, 2, any))+1) %>% 
    write_rds("data/input_processed/lee_species.Rds")
  
  lee_sgb_names %>% 
    write_tsv("data/input_processed/lee_sgb_names.tsv")
  
  lee_metadata %>% 
    filter(Bases >= 1e9) %>% 
    select(sample_id, age_range, PHENOTYPE, `Sample Name`, sex, Bases) %>% 
    write_tsv("data/input_processed/lee_metadata.tsv")
  
  
  ## Thomas et al.
  thomas_mp4 <- read_tsv("data/input_preprocessed/Thomas_2019/profiles.tsv", skip = 1)
  
  thomas_meta <- read_csv("data/input_preprocessed/Thomas_2019/thomas_sra_data.txt", col_types = cols()) %>% 
    rename(sample_id = Run) %>% 
    select(sample_id, `Sample Name`) %>% 
    left_join(read_tsv("data/ext_mg_data/20220322_mg_crc_sample_data.tsv", col_types = cols()) %>% 
                rename(s_id = sample_id,
                       `Sample Name` = subject_id) %>% 
                select(`Sample Name`, study_condition, study_name, disease, disease_subtype, age, gender), by = "Sample Name")
    
  thomas_mp4_sgbs <-
    thomas_mp4 %>% 
    filter(str_detect(clade_name, "t__")) %>% 
    mutate(clade_name = str_remove(clade_name, "^.*t__")) %>% 
    pivot_longer(-clade_name, names_to = "sample_id", values_to = "abundance") %>% 
    group_by(sample_id) %>% 
    mutate(abundance = abundance/sum(abundance)*100) %>% 
    ungroup() %>% 
    pivot_wider(names_from = clade_name, values_from = abundance)
  
  thomas_mp4_sgb_names <- 
    thomas_mp4 %>% 
    filter(str_detect(clade_name, "t__")) %>% 
    select(clade_name) %>% 
    separate(clade_name, into = c("kingdom", "phylum", "clade", "order", "family", "genus", "species", "sgb"), 
             sep = "\\|", 
             remove = FALSE) %>% 
    mutate(across(-clade_name, .fns = function(x) x %>% str_remove("[:alpha:]__")))
  
  thomas_meta %>% 
    write_tsv("data/input_processed/thomas_metadata.tsv")
  thomas_mp4_sgbs %>% 
    write_rds("data/input_processed/thomas_sgbs.Rds")
  thomas_mp4_sgb_names %>% 
    write_tsv("data/input_processed/thomas_sgb_names.tsv")
  
  fit_stool_mp4 <- 
    read_tsv("data/input_preprocessed/fit_stool_it/profiles.tsv", skip = 1, col_types = cols()) %>% 
    pivot_longer(-clade_name) %>% 
    bind_rows(read_tsv("data/input_preprocessed/fit_stool_it/profiles_sub.tsv", skip = 1, col_types = cols()) %>%
                pivot_longer(-clade_name)) %>% 
    pivot_wider(names_from = "name", values_from = "value", values_fill = 0)
  
  fit_stool_meta <- 
    fit_stool_mp4 %>% 
    select(-1) %>% 
    names() %>% 
    enframe(value = "sample_id") %>% 
    select(-1) %>% 
    mutate(sample_type = str_extract(sample_id, "F(E|IT)?") %>% factor(levels = c("FIT", "FE"), labels = c("FIT", "Norgen"))) %>% 
    mutate(id = str_extract(sample_id, "^IT-[:digit:]*")) %>% 
    mutate(subs = case_when(str_detect(sample_id, "(E|IT)$") ~ "original",
                            str_detect(sample_id, "A$") ~ "subs_20e6",
                            str_detect(sample_id, "B$") ~ "subs_15e6",
                            str_detect(sample_id, "C$") ~ "subs_10e6",
                            str_detect(sample_id, "D$") ~ "subs_5e6") %>% 
             factor(levels = c("original", "subs_20e6", "subs_15e6", "subs_10e6", "subs_5e6")))
  
  
  fit_stool_mp4_sgbs <-
    fit_stool_mp4 %>% 
    filter(str_detect(clade_name, "t__")) %>% 
    mutate(clade_name = str_remove(clade_name, "^.*t__")) %>% 
    pivot_longer(-clade_name, names_to = "sample_id", values_to = "abundance") %>% 
    group_by(sample_id) %>% 
    mutate(abundance = abundance/sum(abundance)*100) %>% 
    ungroup() %>% 
    pivot_wider(names_from = clade_name, values_from = abundance)
  
  fit_stool_mp4_sgb_names <- 
    fit_stool_mp4 %>% 
    filter(str_detect(clade_name, "t__")) %>% 
    select(clade_name) %>% 
    separate(clade_name, into = c("kingdom", "phylum", "clade", "order", "family", "genus", "species", "sgb"), 
             sep = "\\|", 
             remove = FALSE) %>% 
    mutate(across(-clade_name, .fns = function(x) x %>% str_remove("[:alpha:]__")))
  
  fit_stool_meta %>% 
    write_tsv("data/input_processed/fit_stool_metadata.tsv")
  fit_stool_mp4_sgbs %>% 
    write_rds("data/input_processed/fit_stool_sgbs.Rds")
  fit_stool_mp4_sgb_names %>% 
    write_tsv("data/input_processed/fit_stool_sgb_names.tsv")
  
}