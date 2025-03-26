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
df_full = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens_contigs/pathogen_species_counts.RDS"))
df_percent = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens_contigs/pathogen_species_percentages.RDS"))
df_full <- data.frame(name = row.names(df_full), df_full)
df_percent_dt <- data.frame(name = row.names(df_percent), df_percent)
setnames(df_percent_dt, old = c('name'), new = c('Species'))

df_percent  <- melt(setDT(df_percent_dt), id.vars = c("Species"), variable.name = "Sample", value.name = "percent")
df_percent_threshold = df_percent[df_percent$percent > 0.001,]

## LONG
rownames(df_full) = df_full$name
setnames(df_full, old = c('name'), new = c('tax_id'))
df_cleaned <- melt(setDT(df_full), id.vars = c("tax_id"), variable.name = "sample_name", value.name = "nt_count")
df_cleaned = df_cleaned[df_cleaned$tax_id %in% df_percent_threshold$Species,]
df_cleaned = merge(df_cleaned, conversion_id, by = "sample_name")

# format_data_diversity_analysis - Species
dat_wide_tax = df_cleaned[, c("tax_id","sample_id", "nt_count")]
dat_wide_tax = dcast(dat_wide_tax, sample_id ~ tax_id, value.var = "nt_count")
dat_wide_tax[is.na(dat_wide_tax)] = 0
dat_wide_tax %<>% as.data.frame()
rownames(dat_wide_tax) = dat_wide_tax$sample_id
dat_wide_tax$sample_id = NULL


# ALPHA DIVERSITY Species
data_richness <- estimateR(dat_wide_tax)
data_evenness <- diversity(dat_wide_tax) / log(specnumber(dat_wide_tax))    
data_shannon <- diversity(dat_wide_tax, index = "shannon") 
data_simpson <- diversity(dat_wide_tax, index = "simpson") 
data_invsimpson <- diversity(dat_wide_tax, index = "invsimpson") 
data_alphadiv_tax <- cbind(t(data_richness),
                           data_shannon,
                           data_simpson,
                           data_invsimpson,
                           data_evenness) %>% as.data.frame()
rm(data_richness, data_evenness, data_shannon, data_invsimpson,data_simpson)  

data_alphadiv_tax$sample_id = rownames(data_alphadiv_tax)
data_alphadiv_tax$sample_type = ifelse(grepl("C", data_alphadiv_tax$sample_id), "Cow dung","Soil floor")
data_alphadiv_tax = merge(data_alphadiv_tax, conversion_id, by = "sample_id")

data_alphadiv_tax <- data_alphadiv_tax[order(data_alphadiv_tax$sample_type, data_alphadiv_tax$hhid),]
data_alphadiv_tax = as.data.table(data_alphadiv_tax)

p1_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$S.obs, y = data_alphadiv_tax[sample_type == "Cow dung",]$S.obs, paired = TRUE, alternative = "two.sided")
P1 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=S.obs)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Richness', x= ' ', y= '', tag = "a)") +
  geom_point() + labs(fill = "Sample Type")  + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p1_stat$p.value, 3)), size = 3)

p2_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$data_evenness, y = data_alphadiv_tax[sample_type == "Cow dung",]$data_evenness, paired = TRUE, alternative = "two.sided")
P2 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=data_evenness)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Evenness', x= ' ', y= '', tag = "b)") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p2_stat$p.value,3)), size = 3)

p3_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$data_shannon, y = data_alphadiv_tax[sample_type == "Cow dung",]$data_shannon, paired = TRUE, alternative = "two.sided")
P3 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=data_shannon)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Shannon', x= ' ', y= '', tag = "c)") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p3_stat$p.value,3)), size = 3)

p4_stat = wilcox.test(x = data_alphadiv_tax[sample_type == "Soil floor",]$data_simpson, y = data_alphadiv_tax[sample_type == "Cow dung",]$data_simpson, paired = TRUE, alternative = "two.sided")
P4 <- ggplot(data_alphadiv_tax, aes(x=sample_type, y=data_simpson)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Simpson', x= ' ', y= '', tag = "d)") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p4_stat$p.value, 2)), size = 3)

aggregate(data_alphadiv_tax$S.obs, list(data_alphadiv_tax$sample_type), FUN=median); p1_stat
aggregate(data_alphadiv_tax$S.chao1, list(data_alphadiv_tax$sample_type), FUN=median); p2_stat
aggregate(data_alphadiv_tax$data_evenness, list(data_alphadiv_tax$sample_type), FUN=median); p3_stat
aggregate(data_alphadiv_tax$data_shannon, list(data_alphadiv_tax$sample_type), FUN=median); p4_stat

combined_species_fig = (P1 | P2) / (P3 | P4) & theme(legend.position = "bottom") 
combined_species_fig = combined_species_fig + plot_layout(guides = "collect")
combined_species_fig

ggsave(
  combined_species_fig,
  filename = paste0(figure_path, "fig-pathogen-alpha-species-bracken", ".pdf"),
  height = 6,
  width = 6.5
)

