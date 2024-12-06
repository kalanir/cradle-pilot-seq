#######################################
# CRADLE trial

# Heatmap of pathogens by sample 
# 2023 environmental pilot
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

species_percent = species_percent %>% as.data.frame()%>% 
  rownames_to_column(var = "species") %>% 
  pivot_longer(-species, names_to = "sample_name", values_to = "percent") %>% 
  mutate(sample_type = ifelse(str_starts(sample_name, "C"), "Cow", "Soil"))

# Fix sample_id mismatch
corrected_ids <- read_excel(paste0(data_dir, "Bracken_processed_data/id_master_list.xlsx")) %>% select(sample_name, sample_id)

species_percent <- left_join(species_percent, corrected_ids, by = "sample_name") %>% select(sample_id, sample_type, species, percent)

# General data visuals & threshold testing -------------------------------------
species_mean_percent <- species_percent %>% 
  group_by(species) %>%
  dplyr::summarize(mean_percent = mean(percent))

# Create histogram with bins of 0.5 increments
hist(species_mean_percent$mean_percent, 
     breaks = seq(floor(min(species_mean_percent$mean_percent)), 
                  ceiling(max(species_mean_percent$mean_percent)), 
                  by = 0.5), 
     main = "Histogram of Species Mean Percentages", 
     xlab = "Mean Percent", 
     ylab = "Frequency")

species_mean_percent <- species_mean_percent %>% 
  filter(mean_percent > 0.5) %>%
  pull(species)


# Make heatmap ----------------------------------------------------------------
heatmap_matrix = species_percent %>%
  filter(species %in% species_mean_percent) %>% 
  pivot_wider(values_from = percent, names_from = species, values_fill = 0, id_cols = sample_id)

## Define the correct sample order
ordered_sample_ids <- c(paste0("C", 1:10), paste0("S", 1:10))

## Set the sample order in heatmap_matrix
heatmap_matrix = heatmap_matrix %>% 
  mutate(sample_id = factor(sample_id, levels = ordered_sample_ids)) %>% 
  arrange(sample_id)

rownames = heatmap_matrix$sample_id
heatmap_matrix = heatmap_matrix %>% dplyr::select(-sample_id)
rownames(heatmap_matrix) = rownames
heatmap_matrix = as.matrix(heatmap_matrix, nrow = 20, ncol = 44) %>% t()
heatmap_matrix = heatmap_matrix + 1

## Log transform data & define a log scale for visualization
log_heatmap_matrix = log10(heatmap_matrix)

# Define the range for the breaks
min_val <- min(log_heatmap_matrix[log_heatmap_matrix > 0], na.rm = TRUE)
max_val <- max(log_heatmap_matrix, na.rm = TRUE)

# Create breaks on a log10 scale
n_breaks <- 7  # Number of breaks (adjust as needed)
log_scale <- seq(min_val, max_val, length.out = n_breaks)
log_scale_norm_transform <- 10^log_scale
log_scale_norm_transform <- round(log_scale_norm_transform, 1)

## Create the color function using log scale
col_fun = colorRamp2(c(0, log_scale[1]-1e-6, log_scale[1:7]), 
                     c("#e8e8e8", "#e8e8e8", "#FED976", "#FEB24C", "#FD8D3C", 
                       "#FC4E2A", "#E31A1C", "#B10026", "#730000"))

pdf(paste0(figure_path, "fig-bracken-pathogen-heatmap.pdf"), width = 6, height = 6.5)

## Run the heatmap
Heatmap(log_heatmap_matrix,
        show_row_names = TRUE,
        show_row_dend = TRUE,
        show_column_dend = FALSE,
        column_split = list(Sample_Slice = cut(1:ncol(heatmap_matrix), breaks = c(0, 10, 20), labels = FALSE)),
        top_annotation = HeatmapAnnotation(
          `sample_type` = anno_block(gp = gpar(fill = c("grey", "#8a4811")),
                                     height = unit(0.5, "cm"),
                                     labels = c("Cow", "Soil"),
                                     labels_gp = gpar(col = "white", fontsize = 9))),
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        column_title = "Sample Type",
        column_title_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 9),
        col = col_fun,
        heatmap_legend_param = list(title = "Relative\nAbundance (%)",
                                    labels = c("Not detected ", log_scale_norm_transform),
                                    at = c(0, log_scale[1:7]),
                                    title_gp = gpar(fontsize = 8),
                                    labels_gp = gpar(fontsize = 8),
                                    legend_height = unit(5, "cm"),
                                    grid_width = unit(0.3, "cm"),
                                    break_dist = c(rep(1, length(log_scale)))
                                    ))

dev.off()

