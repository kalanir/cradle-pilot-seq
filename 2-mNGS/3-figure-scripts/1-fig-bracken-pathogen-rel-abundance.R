#######################################
# CRADLE trial
# 202 environmental pilot

# Sequencing analysis
# Relative abundance plot
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))
library(RColorBrewer)
library(Polychrome)
#remotes::install_github("KarstensLab/microshades")
library(microshades)
library(readxl)

# Read in processed metagenomics data
species_percent = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens_contigs/pathogen_species_percentages.RDS"))
species_percent = species_percent %>% as.data.frame()%>% 
  rownames_to_column(var = "species") %>% 
  pivot_longer(-species, names_to = "sample_name", values_to = "percent") %>% 
  mutate(sample_type = ifelse(str_starts(sample_name, "C"), "Cow", "Soil"))

# Fix sample_id mismatch
corrected_ids <- read_excel(paste0(data_dir, "Bracken_processed_data/id_master_list.xlsx")) %>% select(sample_name, sample_id)

species_percent <- left_join(species_percent, corrected_ids, by = "sample_name") %>% select(sample_id, sample_type, species, percent)

species_percent = species_percent %>% mutate(genus_name = word(species, 1))

# species relative abundance plot ------------------------------------------
# identify 30 most abundant species in both sample types
top_30_species <- species_percent %>% 
  group_by(species) %>% 
  dplyr::summarize(mean_percent = mean(percent)) %>% 
  arrange(desc(mean_percent)) %>% 
  head(30) %>% 
  pull(species)

species_plot_df <- species_percent %>% 
  mutate(species = ifelse(species %in% top_30_species, species, "Other")) %>% 
  group_by(sample_type, sample_id, species) %>% 
  summarise(percent = sum(percent)) %>% 
  mutate(species = factor(species, levels = c(sort(top_30_species), "Other"))) %>%
  mutate(sample_type = ifelse(sample_type == "cow", "Cow Dung", "Soil Floor")) %>% 
  mutate(sample_id = factor(sample_id, levels = c(
    "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"))) %>% 
  mutate(genus_name = word(species, 1))

# define color palette  ------------------------------------------
genus_names = species_plot_df %>% group_by(genus_name, species) %>% dplyr::summarize(count = n()) %>% 
  pull(genus_name) %>% table() %>% as.data.frame()
colnames(genus_names) = c("genus_name", "Freq")

genus_names_multiple = genus_names %>% filter(Freq > 1) %>% pull(genus_name) %>% as.vector()
multiple_color_df = species_plot_df %>% filter(genus_name %in% genus_names_multiple) %>% 
  group_by(genus_name, species) %>% dplyr::summarize(count = n())

## manually assign hex values from microshades palette for those that have multiple species based on counts in genus_names_multiple
microshades_values <-c(microshades_palette("micro_green", n = 3, lightest = F), 
                       microshades_palette("micro_blue", n = 2, lightest = F), 
                       microshades_palette("micro_purple", n = 2, lightest = F),
                       microshades_palette("micro_brown", n = 2, lightest = F),
                       microshades_palette("micro_orange", n = 3, lightest = F),
                       microshades_palette("micro_gray", n = 2, lightest = F))
multiple_color_df = cbind(multiple_color_df, color = microshades_values)

genus_names_single = genus_names %>% filter(Freq == 1) %>% filter(!genus_name == "Other") %>% 
  pull(genus_name) %>% as.vector()
single_color_df = species_plot_df %>% filter(genus_name %in% genus_names_single) %>% 
  group_by(species) %>% 
  dplyr::summarize(count = n())

## create palette of colors sufficiently different from microshades values for the other species
large_palette <- createPalette(N = 100, seedcolors = c("#FFFFFF", "#000000"))

### filter out colors that are too close to the microshades colors
color_distance <- function(color1, color2) {
  sum((col2rgb(color1) - col2rgb(color2))^2)
}

filtered_palette <- large_palette[!sapply(large_palette, function(color) {
  any(sapply(microshades_values, function(ex_color) {
    color_distance(color, ex_color) < 100  # Adjust threshold as needed
  }))
})]
single_colors <- filtered_palette[1:length(genus_names_single)]
swatch(single_colors)

single_color_df = cbind(single_color_df, color = single_colors)

color_df = rbind(multiple_color_df, single_color_df)
other_row = data.frame(species = "Other", Freq = 1, color = "#FFFFFF")
color_df = rbind(color_df, other_row) %>% as.data.frame() %>% select(species, color)

# create final plot
species_plot_df2 = species_plot_df %>% 
  mutate(sample_id = factor(sample_id, levels = c(
    "C1", "S1", "C2", "S2", "C3", "S3", "C4", "S4", "C5", "S5",
    "C6", "S6", "C7", "S7", "C8", "S8", "C9", "S9", "C10", "S10"))) %>% 
  mutate(sample_id_no = as.numeric(gsub("[^0-9]", "", sample_id)))

species_plot_df2 = left_join(species_plot_df2, color_df, by = "species")
colors = setNames(species_plot_df2$color, species_plot_df2$species)

## so that "other" appears last
genus_order = species_plot_df2 %>% filter(!species == 'Other') %>% pull(species) %>% 
  as.vector() %>% unique() %>% sort()
species_plot_df2 = species_plot_df2 %>% mutate(species = factor(species, levels = c(genus_order, "Other")))

plot_species <- ggplot(species_plot_df2, aes(x = sample_id, y = percent, fill = sample_type, 
                                             group = sample_id_no)) + 
  geom_bar(aes(fill = species), stat = "identity", position="fill", col="black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_manual(values = colors) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 25, 50, 75, 100)) + 
  labs(x = "Sample ID", y = "Relative abundance (%)", fill = "Taxon") + 
  facet_grid(~sample_id_no, scales = "free") +
  theme_minimal()+
  ggtitle("b)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color ="black"),
        axis.title.y = element_text(size=14, color="black"),
        plot.title = element_text(size=18)) 

# Make genus plot using only genera from top 30 species ------------------------
genus_plot_df <- species_plot_df2 %>%
  group_by(sample_type, sample_id, genus_name) %>%
  summarise(percent = sum(percent),
            color = first(color)) %>%
  mutate(sample_type = ifelse(sample_type == "cow", "Cow Dung", "Soil Floor"),
         sample_id = factor(sample_id, levels = c(
           "C1", "S1", "C2", "S2", "C3", "S3", "C4", "S4", "C5", "S5",
           "C6", "S6", "C7", "S7", "C8", "S8", "C9", "S9", "C10", "S10"
         )))

genus_plot_df2 = genus_plot_df %>%
  mutate(sample_id_no = as.numeric(gsub("[^0-9]", "", sample_id)))

genus_colors = setNames(genus_plot_df2$color, genus_plot_df2$genus_name)

## so that "other" appears last
genus_order = genus_plot_df2 %>% filter(!genus_name == 'Other') %>% pull(genus_name) %>% 
  as.vector() %>% unique() %>% sort()
genus_plot_df2 = genus_plot_df2 %>% 
  mutate(genus_name = factor(genus_name, levels = c(genus_order, "Other")))

plot_genus <- ggplot(genus_plot_df2, aes(x = sample_id, y = percent, fill = sample_type, 
                                         group = sample_id_no)) + 
  geom_bar(aes(fill = genus_name), stat = "identity", position="fill", col="black") + 
  facet_grid(~sample_id_no, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_manual(values = genus_colors) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 25, 50, 75, 100)) + 
  labs(x = "Sample ID", y = "Relative abundance (%)", fill = "Genus") + 
  theme_minimal() +
  ggtitle("a)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color ="black"),
        axis.title.y = element_text(size=14, color="black"),
        plot.title = element_text(size=18))


plot_both <- grid.arrange(plot_genus, plot_species, ncol=1)

ggsave(plot_both, filename = paste0(figure_path, "fig-bracken-rel-abundance.pdf"), 
       width = 10.2, height = 12.5, dpi = 300)


