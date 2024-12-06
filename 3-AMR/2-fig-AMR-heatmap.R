#######################################
# CRADLE trial

# Visualization of processed sequencing AMR 
# data from 2023 environmental pilot
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))
library(RColorBrewer)
library(circlize)
library(leaflet)
library(ComplexHeatmap)
library(ggh4x)

# read in data
AMR_gene_data = readRDS(paste0(data_dir, 
                              "/CZ_processed_data/AMR/AMR_data_filtered.RDS"))

# gene heatmap ------------------------------------------------------------------
AMR_gene_data = AMR_gene_data %>% 
  mutate(gene = as.character(gene)) %>% 
  mutate(gene = ifelse(gene=="Nocardia farcinica rox","Rox-nf",gene)) %>% 
  mutate(gene = ifelse(gene=="AAC(6')-Ie-APH(2'')-Ia bifunctional protein","AAC(6')-Ie-APH(2'')-Ia",gene)) %>% 
  mutate(gene = ifelse(gene=="vanR gene in vanO cluster", "vanRO", gene)) %>% 
  mutate(gene = paste0(tolower(substr(gene, 1, 1)), substr(gene, 2, nchar(gene)))) %>% 
  mutate(gene = ifelse(gene=="aAC(6')-Ie-APH(2'')-Ia","AAC(6')-Ie-APH(2'')-Ia",gene)) %>% 
  mutate(cutoff_numeric = case_when(
    cutoff == "Perfect" ~ 1,
    cutoff == "Strict" ~ 2,
    cutoff == "Nudged" | cutoff=="Nudged;Strict" | cutoff=="Strict;Nudged" ~ 3,
    TRUE ~ 0
  )) %>% 
  mutate(cutoff_cat = case_when(
    cutoff == "Perfect" ~ "Perfect",
    cutoff == "Strict" | cutoff=="Strict;Nudged" | cutoff=="Nudged;Strict" ~ "Strict",
    cutoff == "Nudged"  ~ "Nudged",
    contigs == 0 ~ "Read"
  ))  %>% 
  mutate(cutoff_cat = factor(cutoff_cat, levels = c("Perfect", "Strict","Nudged", "Read"))) %>% 
  mutate(gene = fct_rev(gene)) %>% 
  mutate(drug_class_cat = case_when(
    drug_class == "aminoglycoside antibiotic" ~ "Aminoglycoside",
    drug_class == "cephamycin" ~ "Cephamycin",
    drug_class == "lincosamide antibiotic" ~ "Lincosamide",
    drug_class == "tetracycline antibiotic" ~ "Tetracycline",
    drug_class == "sulfonamide antibiotic" ~ "Sulfonamide",
    drug_class == "rifamycin antibiotic" ~ "Rifamycin",
    drug_class == "mupirocin-like antibiotic" ~ "Mupirocin",
    drug_class == "diaminopyrimidine antibiotic" ~ "Diaminopyrimidine",
    TRUE ~ "Multiple antibiotics"
  )) %>% 
  mutate(sample_id = factor(sample_id, levels = c(
    "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"
  ))) %>% 
  mutate(sample_type = ifelse(sample_type=="cow", "Cow dung", "Soil floor")) %>% 

  mutate(gene = as.factor(gene)) %>% 
  mutate(gene = fct_rev(gene))

pdf(paste0(figure_path, "fig-AMR-heatmap.pdf"), width=6.5, height=8)
ggplot(AMR_gene_data, aes(x = sample_id, y = gene)) + 
  geom_tile(aes(fill=cutoff_cat), col="black") + 
  scale_fill_manual("Highest\nalignment\nconfidence", values = c("#154696","#6c8ec4","#c6d9f7","#29ab36")) + 
  facet_grid2(drug_class_cat ~ sample_type, scales="free", space="free_y",
              strip = strip_themed(
                background_x = elem_list_rect(fill = c("grey", "#8a4811"))
              )) +
  xlab("Sample") + ylab("Antibiotic resistance gene") + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    strip.text.y = element_text(angle=0),
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text.x = element_text(color="white"),
    strip.background = element_rect(color="black")
  )
dev.off()



# description in text
AMR_gene_data %>% group_by(sample_type, drug_class, gene) %>% 
  summarize(n = n(), .groups = "drop") %>% 
  filter(n > 4) %>% 
  arrange(sample_type, desc(n)) %>% 
  filter(sample_type=="Cow dung") %>% 
  select(drug_class) %>% 
  distinct()
  
AMR_gene_data %>% group_by(sample_type, drug_class, gene) %>% 
  summarize(n = n(), .groups = "drop") %>% 
  filter(n > 4) %>% 
  arrange(sample_type, desc(n)) %>% 
  filter(sample_type=="Soil floor") %>% 
  select(drug_class) %>% 
  distinct()


