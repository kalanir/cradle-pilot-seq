#######################################
# CRADLE trial

# Analysis of sequencing data from 
# 2023 environmental pilot

# DRAFT 
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

library(data.table)
library(magrittr)
library(vegan)
library(patchwork)
library(readxl)
theme_set(theme_bw())

conversion_id <- read_excel(paste0(data_dir, "Bracken_processed_data/id_master_list.xlsx")) %>% select(sample_name, sample_id)
conversion_id$hhid = gsub("S|C", "", conversion_id$sample_id)

# set order/pairs and color
order_type = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10",
               "S1","S2","S3","S4","S5","S6","S7","S8","S9","S10")
order_paired = c("C1","S1","C2","S2","C3","S3","C4","S4","C5","S5","C6","S6","C7","S7","C8","S8","C9","S9","C10","S10")
sample_type_color = c("brown","grey")
names(sample_type_color) = c("Soil floor", "Cow dung")



# Analysis of pathogen data --------------------------------------------
# Read in processed data
taxon_full = readRDS(paste0(data_dir, "Bracken_processed_data/all/species_counts.RDS"))
taxon_percent = readRDS(paste0(data_dir, "Bracken_processed_data/all/species_percentages.RDS"))

taxon_percent_dt = as.data.table(taxon_percent)
taxon_percent_dt$Species = taxon_percent_dt$name
taxon_percent_dt = taxon_percent_dt[,-c(1:2)]
species_percent  <- melt(setDT(taxon_percent_dt), id.vars = c("Species"), variable.name = "Sample")
colnames(species_percent)[3] = "percent"

species_percent_threshold = species_percent[species_percent$percent > 0.001,]
length(unique(species_percent_threshold$Species))

## LONG
rownames(taxon_full) = taxon_full$name
taxon_full = taxon_full[,-c(1:2)]
row_names = rownames(taxon_full)
taxon_full = as.data.frame(taxon_full)
taxon_full$tax_id = row_names
taxon_cleaned <- melt(setDT(taxon_full), id.vars = c("tax_id"), variable.name = "sample_name")
colnames(taxon_cleaned)[3] = "nt_count"
taxon_cleaned = taxon_cleaned[taxon_cleaned$tax_id %in% species_percent_threshold$Species,]
taxon_cleaned = merge(taxon_cleaned, conversion_id, by = "sample_name")

# format_data_diversity_analysis - Species
## utilizing nt_count (raw data)
taxon_cleaned$value = taxon_cleaned$nt_count
dat_count_organism = setDT(taxon_cleaned)[, count_tax := uniqueN(tax_id), by = sample_id] 
dat_count_organism_lim = dat_count_organism[, c("sample_id", "count_tax")] %>% unique()

dat_wide_tax = taxon_cleaned[, c("tax_id","sample_id", "value")]
dat_wide_tax = dcast(dat_wide_tax, sample_id ~ tax_id, value.var = "value")
dat_wide_tax[is.na(dat_wide_tax)] = 0
dat_wide_tax %<>% as.data.frame()
rownames(dat_wide_tax) = dat_wide_tax$sample_id
dat_wide_tax$sample_id = NULL


##### GENUS #######
taxon_full = readRDS(paste0(data_dir, "Bracken_processed_data/all/genus_counts.RDS"))
taxon_percent = readRDS(paste0(data_dir, "Bracken_processed_data/all/genus_percentages.RDS"))

taxon_percent_dt = as.data.table(taxon_percent)
taxon_percent_dt = as.data.table(taxon_percent)
taxon_percent_dt$Genus = taxon_percent_dt$name
taxon_percent_dt = taxon_percent_dt[,-c(1:2)]
genus_percent <- melt(setDT(taxon_percent_dt), id.vars = c("Genus"), variable.name = "Sample")
colnames(genus_percent)[3] = "percent"

genus_percent_threshold = genus_percent[genus_percent$percent > 0.001,]
length(unique(genus_percent_threshold$Genus))

## LONG
rownames(taxon_full) = taxon_full$name
taxon_full = taxon_full[,-c(1:2)]
row_names = rownames(taxon_full)
taxon_full = as.data.table(taxon_full)
taxon_full$Genus = row_names
taxon_cleaned <- melt(setDT(taxon_full), id.vars = c("Genus"), variable.name = "sample_name")
colnames(taxon_cleaned)[3] = "nt_count"
taxon_cleaned = taxon_cleaned[taxon_cleaned$Genus %in% genus_percent_threshold$Genus,]
taxon_cleaned = merge(taxon_cleaned, conversion_id, by = "sample_name")

# format_data_diversity_analysis - Genus
## utilizing nt_count (raw data)
taxon_cleaned$value = taxon_cleaned$nt_count
dat_count_organism = setDT(taxon_cleaned)[, count_tax := uniqueN(Genus), by = sample_id] 
dat_count_organism_lim = dat_count_organism[, c("sample_id", "count_tax")] %>% unique()

dat_wide_genus = taxon_cleaned[, c("Genus","sample_id", "value")]
dat_wide_genus = dcast(dat_wide_genus, sample_id ~ Genus, value.var = "value")
dat_wide_genus[is.na(dat_wide_genus)] = 0
dat_wide_genus %<>% as.data.frame()
rownames(dat_wide_genus) = dat_wide_genus$sample_id
dat_wide_genus$sample_id = NULL





# ALPHA DIVERSITY Species
data_richness <- estimateR(dat_wide_tax)
data_richness
data_evenness <- diversity(dat_wide_tax) / log(specnumber(dat_wide_tax))    
data_shannon <- diversity(dat_wide_tax, index = "shannon") 
data_simpson <- diversity(dat_wide_tax, index = "simpson") 
data_invsimpson <- diversity(dat_wide_tax, index = "invsimpson") 
data_alphadiv_tax <- cbind(t(data_richness), data_shannon, data_simpson, data_invsimpson, data_evenness) %>% as.data.frame() # combine all indices in one data table

rm(data_richness, data_evenness, data_shannon, data_invsimpson,data_simpson)            

data_alphadiv_tax$sample_id = rownames(data_alphadiv_tax)
data_alphadiv_tax$sample_type = ifelse(grepl("C", data_alphadiv_tax$sample_id), "Cow dung","Soil floor")
data_alphadiv_tax = merge(data_alphadiv_tax, conversion_id, by = "sample_id")

data_alphadiv_tax <- data_alphadiv_tax[order(data_alphadiv_tax$sample_type, data_alphadiv_tax$hhid),]
data_alphadiv_tax = as.data.table(data_alphadiv_tax)

p1_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$S.obs, y = data_alphadiv_tax[sample_type == "Cow dung",]$S.obs, paired = TRUE, alternative = "two.sided")
P1 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=S.obs)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Richness', x= ' ', y= '', tag = "A") + ylim(925, 9000) +
  geom_point() + labs(fill = "Sample Type")  + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p1_stat$p.value, 3)), size = 3)
P1
p2_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$S.chao1, y = data_alphadiv_tax[sample_type == "Cow dung",]$S.chao1, paired = TRUE, alternative = "two.sided")
P2 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=S.chao1)) + ylim(925, 9000)+
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p2_stat$p.value, 3)), size = 3)

p3_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$data_evenness, y = data_alphadiv_tax[sample_type == "Cow dung",]$data_evenness, paired = TRUE, alternative = "two.sided")
P3 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=data_evenness)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Evenness', x= ' ', y= '', tag = "C") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p3_stat$p.value,2)), size = 3)

p4_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$data_shannon, y = data_alphadiv_tax[sample_type == "Cow dung",]$data_shannon, paired = TRUE, alternative = "two.sided")
P4 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=data_shannon)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p4_stat$p.value,2)), size = 3)

p5_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$data_simpson, y = data_alphadiv_tax[sample_type == "Cow dung",]$data_simpson, paired = TRUE, alternative = "two.sided")
P5 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=data_simpson)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Simpson', x= ' ', y= '', tag = "E") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p5_stat$p.value, 2)), size = 3)

p6_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$data_invsimpson, y = data_alphadiv_tax[sample_type == "Cow dung",]$data_invsimpson, paired = TRUE, alternative = "two.sided")
P6 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=data_invsimpson)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Inv Simpson', x= ' ', y= '', tag = "F") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ", round(p6_stat$p.value,2)), size = 3)

aggregate(data_alphadiv_tax$S.obs, list(data_alphadiv_tax$sample_type), FUN=median); p1_stat
aggregate(data_alphadiv_tax$S.chao1, list(data_alphadiv_tax$sample_type), FUN=median); p2_stat
aggregate(data_alphadiv_tax$data_evenness, list(data_alphadiv_tax$sample_type), FUN=median); p3_stat
aggregate(data_alphadiv_tax$data_shannon, list(data_alphadiv_tax$sample_type), FUN=median); p4_stat
aggregate(data_alphadiv_tax$data_simpson, list(data_alphadiv_tax$sample_type), FUN=median); p5_stat
aggregate(data_alphadiv_tax$data_invsimpson, list(data_alphadiv_tax$sample_type), FUN=median); p6_stat


combined = (P1 | P2 | P3) / (P4 | P5 | P6)  & theme(legend.position = "bottom") 
combined_alpha_fig = combined + plot_layout(guides = "collect")
combined_alpha_fig
ggsave(combined_alpha_fig, filename = paste0(figure_path, "fig-all-alpha-species-bracken",".png"), height = 5, width = 8)

# ALPHA DIVERSITY GENUS
data_richness <- estimateR(dat_wide_genus)
data_richness
data_evenness <- diversity(dat_wide_genus) / log(specnumber(dat_wide_genus))    
data_shannon <- diversity(dat_wide_genus, index = "shannon") 
data_simpson <- diversity(dat_wide_genus, index = "simpson") 
data_invsimpson <- diversity(dat_wide_genus, index = "invsimpson") 
data_alphadiv_genus <- cbind(t(data_richness), data_shannon, data_simpson, data_invsimpson, data_evenness) %>% as.data.frame() # combine all indices in one data table

rm(data_richness, data_evenness, data_shannon, data_invsimpson,data_simpson)            

data_alphadiv_genus$sample_id = rownames(data_alphadiv_genus)
data_alphadiv_genus$sample_type = ifelse(grepl("C", data_alphadiv_genus$sample), "Cow dung", "Soil floor")

data_alphadiv_genus = merge(data_alphadiv_genus, conversion_id, by = "sample_id")


data_alphadiv_genus <- data_alphadiv_genus[order(data_alphadiv_genus$sample_type, data_alphadiv_genus$hhid),]
data_alphadiv_genus = as.data.table(data_alphadiv_genus)

p1_stat = wilcox.test(x = data_alphadiv_genus[sample_type == "Soil floor",]$S.obs, y = data_alphadiv_genus[sample_type == "Cow dung",]$S.obs, paired = TRUE, alternative = "two.sided")
P1 <- ggplot(data_alphadiv_genus, aes(x=sample_type, y=S.obs)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Richness', x= ' ', y= '', tag = "A") +
  geom_point() + labs(fill = "Sample Type")  + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p1_stat$p.value, 2)), size = 3)

p2_stat = wilcox.test(x = data_alphadiv_genus[sample_type == "Soil floor",]$S.chao1, y = data_alphadiv_genus[sample_type == "Cow dung",]$S.chao1, paired = TRUE, alternative = "two.sided")
P2 <- ggplot(data_alphadiv_genus, aes(x=sample_type, y=S.chao1)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p2_stat$p.value, 2)), size = 3)

p3_stat = wilcox.test(x = data_alphadiv_genus[sample_type == "Soil floor",]$data_evenness, y = data_alphadiv_genus[sample_type == "Cow dung",]$data_evenness, paired = TRUE, alternative = "two.sided")
P3 <- ggplot(data_alphadiv_genus, aes(x=sample_type, y=data_evenness)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Evenness', x= ' ', y= '', tag = "C") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p3_stat$p.value,2)), size = 3)

p4_stat = wilcox.test(x = data_alphadiv_genus[sample_type == "Soil floor",]$data_shannon, y = data_alphadiv_genus[sample_type == "Cow dung",]$data_shannon, paired = TRUE, alternative = "two.sided")
P4 <- ggplot(data_alphadiv_genus, aes(x=sample_type, y=data_shannon)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p4_stat$p.value,3)), size = 3)

p5_stat = wilcox.test(x = data_alphadiv_genus[sample_type == "Soil floor",]$data_simpson, y = data_alphadiv_genus[sample_type == "Cow dung",]$data_simpson, paired = TRUE, alternative = "two.sided")
P5 <- ggplot(data_alphadiv_genus, aes(x=sample_type, y=data_simpson)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Simpson', x= ' ', y= '', tag = "E") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ", round(p6_stat$p.value,3)), size = 3)

p6_stat = wilcox.test(x = data_alphadiv_genus[sample_type == "Soil floor",]$data_invsimpson, y = data_alphadiv_genus[sample_type == "Cow dung",]$data_invsimpson, paired = TRUE, alternative = "two.sided")
P6 <- ggplot(data_alphadiv_genus, aes(x=sample_type, y=data_invsimpson)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Inv Simpson', x= ' ', y= '', tag = "F") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ", round(p6_stat$p.value,3)), size = 3)

combined = (P1 | P2 | P3) / (P4 | P5 | P6) & theme(legend.position = "bottom") 
combined_alpha_fig = combined + plot_layout(guides = "collect")
combined_alpha_fig
ggsave(combined_alpha_fig, filename = paste0(figure_path, "fig-all-alpha-genus-bracken",".png"), height = 5, width = 8)

# END

