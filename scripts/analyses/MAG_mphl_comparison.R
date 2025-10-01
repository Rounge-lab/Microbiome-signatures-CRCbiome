

## Metaphlan 4
mp4_gtdb <- read_tsv("data/input_preprocessed/gtdb_profiles.tsv", skip = 1)

mp4_gtdb_species <-
  mp4_gtdb %>% 
  filter(str_detect(clade_name, "s__")) %>% 
  (function(x) {
    bind_cols(x  %>% select(clade_name), 
              x %>% 
                select(clade_name) %>% 
                separate_wider_delim(cols = clade_name, delim = ";", 
                                     names = c("domain", 
                                               "phylum",
                                               "clade",
                                               "order",
                                               "family",
                                               "genus",
                                               "species")),
              x %>% select(-clade_name)
    )
  }) %>% 
  filter(!species == "s__") %>% 
  select(species, starts_with("S_")) %>% 
  mutate(species = str_remove(species, "s__")) %>% 
  pivot_longer(-species, values_to = "mphl_abundance", names_to = "sample_id")
  

MAGs_species <-
  MAGs_abundance %>% 
  pivot_longer(-sample_id, names_to = "MAG_id", values_to = "MAG_abundance") %>% 
  inner_join(MAG_taxonomy %>% select(MAG_id, species) %>% filter(!is.na(species))) %>% 
  group_by(sample_id, species) %>%
  summarize(MAG_abundance = mean(MAG_abundance),
            MAG_ids = paste0(unique(MAG_id), collapse = ","),
            .groups = "drop")

## What is the overlap?
MAGs_species %>% 
  distinct(species) %>% 
  mutate(in_mphl = species %in% mp4_gtdb_species$species) %>% 
  count(in_mphl)
mp4_gtdb_species %>% 
  distinct(species) %>% 
  mutate(in_mags = species %in% MAGs_species$species) %>% 
  count(in_mags)

## Calculate cosine correlation for each species 
comparison_dat <-
  MAGs_species %>% 
  inner_join(mp4_gtdb_species, by = join_by(species, sample_id)) %>% 
  group_by(species) %>% 
  filter(any(MAG_abundance > 0), any(mphl_abundance > 0)) %>% 
  summarize(cosine_cor = lsa::cosine(MAG_abundance, mphl_abundance)[,1],
            pearson_cor = cor(MAG_abundance, mphl_abundance),
            prev_mags = mean(MAG_abundance > 0),
            prev_mphl = mean(mphl_abundance > 0),
            n_mags = str_count(unique(MAG_ids), "MAG"),
            .groups = "drop") %>%
  arrange(cosine_cor) 

cor.test(log10(comparison_dat$prev_mags+1e-5), log10(comparison_dat$prev_mphl+1e-5))
cor.test(comparison_dat$prev_mags, comparison_dat$prev_mphl, method = "spearman")
comparison_dat %>% 
  skimr::skim()

comparison_dat %>% 
  mutate(high_cosine = cosine_cor > 0.95) %>% 
  count(high_cosine) %>% 
  mutate(frac = n/sum(n)*100)

comparison_plot <-
  comparison_dat %>% 
  ggplot(aes(x = prev_mags, y = prev_mphl, color = cosine_cor)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(x = "MAG prevalence",
       y = "Metaphlan prevalence",
       color = "Cosine similarity\n measure") +
  theme(text = element_text(size = 6),
        legend.key.size = unit(4, units = "mm"))

comparison_plot %>% 
  ggsave2(filename = "results/figures/misc/comparison_metaphlan_MAGs.pdf", plot = ., height = 70, width = 90, units = "mm")
