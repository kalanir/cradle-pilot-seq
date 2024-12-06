#######################################
# CRADLE trial

# Visualization of processed sequencing AMR 
# data from 2023 environmental pilot
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))
library(Polychrome)

# read in data
AMR_filtered = readRDS(paste0(data_dir, 
                              "/CZ_processed_data/AMR/AMR_data_filtered.RDS"))

drug_class_by_sample_type <- AMR_filtered %>%
  group_by(sample_type, hhid, drug_class_cat) %>%
  summarize(num_genes = n(), .groups = "drop") %>% 
  group_by(sample_type, drug_class_cat) %>%
  summarize(num_genes = sum(num_genes)) %>% 
  ungroup() %>% 
  complete(sample_type, drug_class_cat, fill = list(num_genes = 0)) %>%
  mutate(sample_type = ifelse(sample_type == "cow", "Cow Dung", "Soil Floor")) %>% 
  mutate(drug_class_cat = ifelse(drug_class_cat == "Disinfecting Agents And Antiseptics","Disinfecting agents and antiseptics", drug_class_cat))

# Calculate the frequency of each group within each category
df_freq <- drug_class_by_sample_type %>%
  group_by(drug_class_cat) %>%
  summarise(N_reads = sum(num_genes)) %>%
  ungroup() %>% 
  arrange(N_reads)

#----------------------------------------------------------------
# make a bar plot showing drug classes by sample type
# dodge by drug class
#----------------------------------------------------------------

bar_plot <- ggplot(drug_class_by_sample_type, aes(x = factor(drug_class_cat, levels = df_freq$drug_class_cat), 
                                                  y = num_genes, fill = sample_type)) +
  geom_bar(stat = "identity", position = "dodge",width=0.75) +
  theme_minimal() +
  labs(x = "Drug class", y = "Number of AMR genes", fill = "Sample type") +
  scale_fill_manual(values = c("grey", "#8a4811")) +
  scale_y_continuous(breaks = seq(0,60,5),
                     labels = seq(0,60,5)) +
  coord_flip() + 
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_line(linewidth = 0.1),
        panel.grid.major.y = element_blank())
bar_plot

ggsave(paste0(figure_path, "fig-AMR-drug_class_barplot.pdf"), bar_plot, width = 6, height = 5)


#----------------------------------------------------------------
# make a stacked bar plot with the individual drug classes
# as colors and x-axis as sample type and then y axis as number
# of genes
#----------------------------------------------------------------
multi <- AMR_filtered %>%
  filter(drug_class_cat=="Multiple antibiotics") 

# separate individual drugs within drug_class column
multi <- multi %>% 
  separate_rows(drug_class, sep = ";") %>% 
  mutate(drug_class = sub("^\\s+", "", drug_class)) %>% 
  mutate(drug_class = str_to_title(drug_class)) %>% 
  filter(drug_class != "Multiple antibiotics") %>% 
  group_by(sample_type, drug_class) %>% 
  summarise(num_genes = n(), .groups = "drop") %>% 
  mutate(drug_class = fct_reorder(drug_class, num_genes, .fun = sum)) %>% 
  mutate(sample_type = ifelse(sample_type == "cow", "Cow Dung", "Soil Floor")) %>% 
  mutate(drug_class = as.character(drug_class)) %>% 
  mutate(drug_class = ifelse(drug_class == "Disinfecting Agents And Antiseptics","Disinfecting agents and antiseptics", drug_class)) %>% 
  mutate(drug_class = ifelse(drug_class == "Penam", "Penem", drug_class)) %>% 
  mutate(drug_class = gsub(" Antibiotic", "", drug_class)) 
  
mypal <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#FFFF33",  # Yellow
  "#A65628",  # Brown
  "#F781BF",  # Pink
  "#1B9E77",  # Teal
  "#D95F02",  # Dark orange
  "#7570B3",  # Blue-purple
  "#E7298A",  # Magenta
  "#66A61E",  # Lime green
  "#E6AB02",  # Gold
  "#A6761D",  # Dark brown
  "#666666",  # Gray
  "#1F77B4",  # Steel blue
  "#FF7F0E",  # Dark yellow
  "#2CA02C",  # Forest green
  "#D62728",  # Crimson
  "#9467BD",  # Light purple
  "#8C564B",  # Sienna
  "#E377C2",  # Light pink
  "#7F7F7F",  # Medium gray
  "#BCBD22",  # Olive
  "#17BECF",  # Cyan
  "#AA40FC",  # Bright purple
  "#B15928",  # Rust
  "#6B6ECF",  # Periwinkle
  "#D4B9DA",  # Light lavender
  "#AEC7E8",  # Light blue
  "#98DF8A"   # Light green
)

# First, calculate the sum of genes for each drug class
drug_totals <- multi %>%
  group_by(drug_class) %>%
  summarise(total_genes = sum(num_genes)) %>%
  arrange(total_genes)  

multi = multi %>%  mutate(drug_class = factor(drug_class, 
                                    levels = drug_totals$drug_class))

multi_plot <- ggplot(multi, aes(x = sample_type, y = num_genes)) +
  geom_bar(aes(fill = drug_class), stat = "identity", col = "black") +
  theme_minimal() +
  labs(x = "Sample type", y = "Number of AMR genes", fill = "Drug class")+
  scale_fill_manual(values = mypal) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size=11, color="black")) 

ggsave(paste0(figure_path, "fig-AMR-drug_class_stacked_barplot.pdf"), multi_plot,
       width = 8, height = 5)


