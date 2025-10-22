
sample_data <- 
  read_rds("data/synthetic/sample_data.Rds")

screening_data <- read_rds("data/synthetic/screening_data.Rds")

train_test <- read_csv("data/synthetic/train_test.csv", col_types = cols())

mphlan_abundance <- read_rds("data/synthetic/mp4_sgbs.Rds")
MAGs_abundance <- read_rds("data/synthetic/MAG_abundance.Rds")

MAG_taxonomy <- read_tsv("data/synthetic/MAG_taxonomy.tsv", col_types = cols())
mp4_taxonomy <- read_tsv("data/synthetic/mp4_sgb_names.tsv", col_types = cols())

ko_hum <- read_rds("data/synthetic/kegg.Rds")

ko_annotations <- read_tsv("data/synthetic/ko_annotations.tsv", col_types = cols())

variables <- read_rds("data/synthetic/metadata_variables.Rds")
meta_dat <- read_rds("data/synthetic/metadata_selected_variables.Rds")
meta_dat_cat <- read_rds("data/synthetic/metadata_selected_variables_cat.Rds")

lesion_locs <- read_tsv("data/synthetic/lesion_locs.tsv", col_types = cols())
