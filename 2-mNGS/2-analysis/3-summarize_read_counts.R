source(here::here("0-config.R"))
library(readr)
library(readxl)

# Fix sample_id mismatch
corrected_ids <- read_excel(paste0(data_dir, "Bracken_processed_data/id_master_list.xlsx")) %>% select(sample_name, sample_id)
sample_readcounts = read_tsv(paste0(data_dir, "Bracken_processed_data/sample_readcounts.tsv")) %>% 
  left_join(corrected_ids, by = c("Sample" = "sample_name")) %>% 
  mutate(Sample = sample_id) %>% 
  select(-sample_id) %>% 
  mutate(sample_group = ifelse(str_detect(Sample, "C"), "Cow", "Soil"),
         Sample = factor(Sample, levels = c(paste0("C", 1:10), paste0("S", 1:10)))) 

readcounts_by_group = sample_readcounts %>% 
  group_by(sample_group) %>% 
  summarize(min_raw_reads = min(raw_reads) / 1000000, 
            max_raw_reads = max(raw_reads) / 1000000, 
            mean_raw_reads = sum(raw_reads) / 1000000, 
            mean_processed_reads = mean(orphan_reads) / 1000000,
            min_processed_reads = min(orphan_reads) / 1000000, 
            max_processed_reads = max(orphan_reads) / 1000000,
            mean_prop_retained = mean(orphan_frac) * 100)

readcounts_by_group

table_a2 = sample_readcounts %>% 
  arrange(Sample) %>% 
  mutate(raw_reads = round(raw_reads / 1000000, 3),
         orphan_reads = round(orphan_reads / 1000000, 3),
         orphan_frac = round(orphan_frac * 100, 2)) %>% 
  select(Sample, 
         `Raw Reads (millions)` = raw_reads, 
         `Processed Reads (millions)` = orphan_reads, 
         `Percent Retained (%)` = orphan_frac)
write_csv(table_a2, here::here("tables", "TableA2_sample_readcounts.csv"))



