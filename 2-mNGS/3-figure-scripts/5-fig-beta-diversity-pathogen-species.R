#######################################
# CRADLE trial

# Analysis of sequencing data from 
# 2025
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

library(data.table)
library(magrittr)
library(vegan)
library(patchwork)
library(paletteer)
library(readxl)
theme_set(theme_bw())

path_cradle_data = paste0(here::here(), "/data/all/bracken/")

# set order/pairs and color
order_type = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10",
               "S1","S2","S3","S4","S5","S6","S7","S8","S9","S10")
order_paired = c("C1","S1","C2","S2","C3","S3","C4","S4","C5","S5","C6","S6","C7","S7","C8","S8","C9","S9","C10","S10")
sample_type_color = c("brown","grey")
names(sample_type_color) = c("Soil floor", "Cow dung")

## conversions
conversion_id <- read_excel(paste0(data_dir, "Bracken_processed_data/id_master_list.xlsx")) %>% select(sample_name, sample_id)
conversion_id$hhid = gsub("S|C", "", conversion_id$sample_id)

# # Analysis of pathogen data --------------------------------------------
metadata = fread(paste0(table_path, "table_A1.csv"))
metadata$hhid %<>% as.character()

metadata = merge(metadata, conversion_id, by = "hhid")
metadata$hhid = factor(metadata$hhid, levels = c(1:10))
metadata$sample_type = ifelse(grepl("C", metadata$sample_id), "Cow dung", "Soil floor")

#### CLEAN INPUT
# Analysis of pathogen data ---------------------------------------------
# Read in processed data
taxon_percent = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens_contigs/pathogen_species_percentages.RDS"))
taxon_percent_dt <- data.frame(name = row.names(taxon_percent), taxon_percent)
setnames(taxon_percent_dt, old = c('name'), new = c('Species'))
species_percent <- melt(setDT(taxon_percent_dt), id.vars = c("Species"), variable.name = "Sample", value.name = "percent")

species_percent_threshold = species_percent[species_percent$percent > 0.001,]
length(unique(species_percent_threshold$Species))

taxon_cleaned = species_percent_threshold
colnames(taxon_cleaned) = c("tax_id","sample_name","percent")

# merge taxon with household data
taxon_cleaned = merge(taxon_cleaned, metadata[,c("sample_name","sample_id")], by = c("sample_name"))

## Beta-diversity using species
dat_wide_tax = taxon_cleaned[, c("tax_id","sample_id", "percent")]
dat_wide_tax = dcast(dat_wide_tax, sample_id ~ tax_id, value.var = "percent")
dat_wide_tax[is.na(dat_wide_tax)] = 0
dat_wide_tax %<>% as.data.frame()
rownames(dat_wide_tax) = dat_wide_tax$sample_id
dat_wide_tax$sample_id = NULL


## set-up metadata
metadata = as.data.frame(metadata)
rownames(metadata) = metadata$sample_id
metadata = metadata[match(rownames(dat_wide_tax), rownames(metadata)),]
all.equal(rownames(metadata), rownames(dat_wide_tax))

# Calculate Bray-Curtis distance among samples and convert the result to a matrix
bray_curtis_dist <- (vegdist(dat_wide_tax, method = "bray"))
dist_bc <- as.matrix(bray_curtis_dist)

# calculate PCOA
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)

# All components could be found here: 
# bray_curtis_pcoa$vectors but look at the first 2
all.equal(rownames(metadata), rownames(bray_curtis_pcoa$vectors))
bray_curtis_pcoa_df <- data.frame(samples = rownames(bray_curtis_pcoa$vectors),pcoa1 = bray_curtis_pcoa$vectors[,1], pcoa2 = bray_curtis_pcoa$vectors[,2], metadata)

## PERMANOVA tests if the centroids, similar to means, of each group are significantly different from each other.
all.equal(rownames(dist_bc), rownames(metadata))
set.seed(4)
p1_stat = adonis2(dist_bc~hhid,data=metadata, permutations=9999, method="bray")
set.seed(4)
p2_stat = adonis2(dist_bc~sample_type,data=metadata, permutations=9999, method="bray")
set.seed(4)
p3_stat = adonis2(dist_bc~dung_floor,data=metadata, permutations=9999, method="bray")


p1_stat = p1_stat$`Pr(>F)`[1]
p2_stat = p2_stat$`Pr(>F)`[1]
p3_stat = p3_stat$`Pr(>F)`[1]
p1_stat; p2_stat; p3_stat

# Create a plot
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = hhid)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA: HHID", color = "Household ID") +
  theme(title = element_text(size = 8.5)) + annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, label = paste0("PERMANOVA = ",p1_stat), size = 3)+ scale_color_manual(values = paletteer_d("ggthemes::Classic_10"))

bray_curtis_plot


bray_curtis_plot2 <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = sample_type)) +
  geom_point()+ scale_color_manual(values = sample_type_color) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA: Sample Type", color = "Sample\nType") +
  theme(title = element_text(size = 8.5)) + annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, label = paste0("PERMANOVA = ",p2_stat), size = 3) 

bray_curtis_plot2


bray_curtis_plot3 <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = dung_floor)) +
  geom_point()+ scale_color_manual(values = c("grey", "#625377FF")) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA: Animal Feces", color =  "Animal Feces\nPresent on\nFloor") +
  theme(title = element_text(size = 8.5))+ annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, label = paste0("PERMANOVA = ",p3_stat), size = 3)

bray_curtis_plot3
combined_species = (bray_curtis_plot | bray_curtis_plot2 | bray_curtis_plot3) & theme(legend.position = "bottom") 
combined_species
ggsave(
  combined_species,
  filename = paste0(figure_path, "fig-pathogen-beta-species-bracken", ".pdf"),
  height = 4.5,
  width = 10
)
