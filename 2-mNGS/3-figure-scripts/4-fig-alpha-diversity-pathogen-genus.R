
#######################################
# CRADLE trial

# Analysis of sequencing data from 
# 2025 environmental pilot

#######################################
rm(list=ls())
# library(devtools)
# install_version("patchwork", version = "1.3.0", repos = "http://cran.us.r-project.org")

source(paste0(here::here(), "/0-config.R"))

library(data.table)
library(magrittr)
library(vegan)
library(patchwork)
library(readxl)
theme_set(theme_bw())

conversion_id <- read_excel(paste0(data_dir, "Bracken_processed_data/id_master_list.xlsx")) %>% select(sample_name, sample_id, hhid)

# set order/pairs and color
order_type = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10",
               "S1","S2","S3","S4","S5","S6","S7","S8","S9","S10")
order_paired = c("C1","S1","C2","S2","C3","S3","C4","S4","C5","S5","C6","S6","C7","S7","C8","S8","C9","S9","C10","S10")
sample_type_color = c("brown","grey")
names(sample_type_color) = c("Soil floor", "Cow dung")

##### GENUS #######
df_full = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens_contigs/pathogen_genus_counts.RDS"))
df_percent = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens_contigs/pathogen_genus_percentages.RDS"))
df_full <- data.frame(name = row.names(df_full), df_full)
df_percent_dt <- data.frame(name = row.names(df_percent), df_percent)
setnames(df_percent_dt, old = c('name'), new = c('Genus'))

df_percent <- melt(setDT(df_percent_dt), id.vars = c("Genus"), variable.name = "Sample", value.name = "percent")
df_percent_threshold = df_percent[df_percent$percent > 0.001,]

## LONG
setnames(df_full, old = c('name'), new = c('Genus'))
df_cleaned <- melt(setDT(df_full), id.vars = c("Genus"), variable.name = "sample_name", value.name = "nt_count")
df_cleaned = df_cleaned[df_cleaned$Genus %in% df_percent_threshold$Genus,]
df_cleaned = merge(df_cleaned, conversion_id, by = "sample_name")

# format_data_diversity_analysis - Genus
dat_wide_genus = df_cleaned[, c("Genus","sample_id", "nt_count")]
dat_wide_genus = dcast(dat_wide_genus, sample_id ~ Genus, value.var = "nt_count")
dat_wide_genus[is.na(dat_wide_genus)] = 0
dat_wide_genus %<>% as.data.frame()
rownames(dat_wide_genus) = dat_wide_genus$sample_id
dat_wide_genus$sample_id = NULL



# ALPHA DIVERSITY GENUS
data_richness <- estimateR(dat_wide_genus)
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
  labs(title= 'Richness', x= ' ', y= '', tag = "a") +
  geom_point() + labs(fill = "Sample Type")  + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p1_stat$p.value, 3)), size = 3)

p2_stat = wilcox.test(x = data_alphadiv_genus[sample_type == "Soil floor",]$data_evenness, y = data_alphadiv_genus[sample_type == "Cow dung",]$data_evenness, paired = TRUE, alternative = "two.sided")
P2 <- ggplot(data_alphadiv_genus, aes(x=sample_type, y=data_evenness)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Evenness', x= ' ', y= '', tag = "b") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p2_stat$p.value,3)), size = 3)

p3_stat = wilcox.test(x = data_alphadiv_genus[sample_type == "Soil floor",]$data_shannon, y = data_alphadiv_genus[sample_type == "Cow dung",]$data_shannon, paired = TRUE, alternative = "two.sided")
P3 <- ggplot(data_alphadiv_genus, aes(x=sample_type, y=data_shannon)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Shannon', x= ' ', y= '', tag = "c") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p3_stat$p.value,3)), size = 3)

p4_stat = wilcox.test(x = data_alphadiv_genus[sample_type == "Soil floor",]$data_simpson, y = data_alphadiv_genus[sample_type == "Cow dung",]$data_simpson, paired = TRUE, alternative = "two.sided")
P4 <- ggplot(data_alphadiv_genus, aes(x=sample_type, y=data_simpson)) +
  geom_boxplot(aes(fill = sample_type)) + scale_fill_manual(values = sample_type_color) + geom_line(aes(group = hhid), color = "#707B7C") +
  labs(title= 'Simpson', x= ' ', y= '', tag = "d") +
  geom_point() + labs(fill = "Sample Type") + annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1.0, label = paste0("p = ",round(p4_stat$p.value, 3)), size = 3)


combined_genus_fig = (P1 | P2) / (P3 | P4) & theme(legend.position = "bottom") 
combined_genus_fig = combined_genus_fig + plot_layout(guides = "collect")
combined_genus_fig


ggsave(
  combined_genus_fig,
  filename = paste0(figure_path, "fig-pathogen-alpha-genus-bracken", ".pdf"),
  height = 6,
  width = 6.5
)

