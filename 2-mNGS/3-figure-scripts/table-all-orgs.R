#######################################
# CRADLE trial

# 2023 environmental pilot
# table of all microbes
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))
library(RColorBrewer)
library(readxl)
library(circlize)
library(leaflet)
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

select <- dplyr::select
summarize <- dplyr::summarize

# Read in processed data ------------------------------------------------------
species_percent = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens/pathogen_species_percentages.RDS"))
rownames(species_percent) = species_percent$name
species_percent = species_percent[,-c(1:2)]

genus_percent = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens/pathogen_genus_percentages.RDS"))

species_percent = species_percent %>% as.data.frame()%>% 
  rownames_to_column(var = "species") %>% 
  pivot_longer(-species, names_to = "sample_name", values_to = "percent") %>% 
  mutate(sample_type = ifelse(str_starts(sample_name, "C"), "Cow", "Soil"))

genus_percent = genus_percent %>% as.data.frame() %>% 
  rename(genus = name) %>%
  dplyr::select(-taxid) %>% 
  pivot_longer(-genus, names_to = "sample_name", values_to = "percent") %>% 
  mutate(sample_type = ifelse(str_starts(sample_name, "C"), "Cow", "Soil"))

# Fix sample_id mismatch
corrected_ids <- read_excel(paste0(data_dir, "Bracken_processed_data/id_master_list.xlsx")) %>% select(sample_name, sample_id)
species_percent <- left_join(species_percent, corrected_ids, by = "sample_name") %>% select(sample_id, sample_type, species, percent)
genus_percent <- left_join(genus_percent, corrected_ids, by = "sample_name") %>% select(sample_id, sample_type, genus, percent)

# Appendix table: species -----------------------------------------------------
table_cow = species_percent %>% filter(percent>0 & sample_type=="Cow") %>% 
  dplyr::select(species) %>% 
  distinct() %>% 
  arrange(species)

table_soil = species_percent %>% filter(percent>0 & sample_type=="Soil") %>% 
  dplyr::select(species) %>% 
  distinct() %>% 
  arrange(species)

write.csv(table_cow, paste0(table_path, "table_cow_species.csv"))
write.csv(table_soil, paste0(table_path, "table_soil_species.csv"))

# Appendix table: genus -------------------------------------------------------
table_cow_genus = genus_percent %>% filter(percent>0 & sample_type=="Cow") %>% 
  dplyr::select(genus) %>% 
  distinct() %>% 
  arrange(genus)

table_soil_genus = genus_percent %>% filter(percent>0 & sample_type=="Soil") %>%
  dplyr::select(genus) %>% 
  distinct() %>% 
  arrange(genus)

write.csv(table_cow_genus, paste0(table_path, "table_cow_genus.csv"))
write.csv(table_soil_genus, paste0(table_path, "table_soil_genus.csv"))
