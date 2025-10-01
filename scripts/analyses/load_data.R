
sample_data <- 
  read_rds("data/input_processed/sample_data.Rds")

screening_data <- read_rds("data/input_processed/screening_data.Rds")

train_test <- read_csv("data/train_test/train_test.csv", col_types = cols())

mphlan_abundance <- read_rds("data/input_processed/mp4_sgbs.Rds")
MAGs_abundance <- read_rds("data/input_processed/MAG_abundance.Rds")

MAG_taxonomy <- read_tsv("data/input_processed/MAG_taxonomy.tsv", col_types = cols())
mp4_taxonomy <- read_tsv("data/input_processed/mp4_sgb_names.tsv", col_types = cols())

ko_hum <- read_rds("data/input_processed/humann_rerun/kegg.Rds")

ko_annotations <- read_tsv("data/input_processed/humann_rerun/ko_annotations.tsv", col_types = cols())

variables <- read_rds("data/input_processed/metadata_variables.Rds")
meta_dat <- read_rds("data/input_processed/metadata_selected_variables.Rds")
meta_dat_cat <- read_rds("data/input_processed/metadata_selected_variables_cat.Rds")

