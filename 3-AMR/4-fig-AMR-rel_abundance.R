#######################################
# CRADLE trial

# Relative abundance of ARGs by drug class
# and by AMR mechanism 

# data from 2023 environmental pilot
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))
library(Polychrome)

# read in data
AMR_filtered = readRDS(paste0(data_dir, 
                              "/CZ_processed_data/AMR/AMR_data_filtered.RDS"))

drug_class_ord <- AMR_filtered %>%
  group_by(drug_class_cat) %>%
  summarize(num_genes = n(), .groups = "drop") %>%
  arrange(desc(num_genes))


# plot of relative abundance of drug classes by sample type ---------------------

drug_class_by_sample_type <- AMR_filtered %>%
  group_by(sample_type, drug_class_cat) %>%
  summarize(num_genes = n(), .groups = "drop") %>% 
  mutate(drug_class_cat = factor(drug_class_cat, levels = 
                                   drug_class_ord$drug_class_cat)) %>% 
  mutate(drug_class_cat = fct_rev(drug_class_cat)) %>% 
  mutate(sample_type = ifelse(sample_type == "cow", "Cow Dung", "Soil Floor"))

qualitative_colors <- c(
  "#4e79a7",  # Blue
  "#f28e2b",  # Orange
  "#59a14f",  # Green
  "#e15759",  # Red
  "#76b7b2",  # Teal
  "#edc948",  # Yellow
  "#b07aa1",  # Purple
  "#9c755f",  # Brown
  "#bab0ac"   # Gray
)

arg_plot <- ggplot(drug_class_by_sample_type, aes(x = sample_type, y = num_genes)) + 
  geom_bar(aes(fill = drug_class_cat), stat = "identity", position="fill", col="black",
           width=0.8) +
  scale_fill_manual(values = qualitative_colors)+
  theme_minimal() +
  scale_y_continuous(labels = c("0","25","50","75","100"),
                     breaks = c(0,.25, .5, .75, 1)) + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black",size=12)) +
  ylab("Relative abundance of ARGs (%)")

ggsave(paste0(figure_path, "fig-AMR-drug_class_rel_abundance.pdf"),
       arg_plot, width = 5, height = 3)


# plot of relative abundance of mechanism by sample type ---------------------
mechanism_ord <- AMR_filtered %>%
  group_by(mechanism) %>%
  summarize(num_genes = n(), .groups = "drop") %>%
  arrange(desc(num_genes))

mechanism_by_sample_type <- AMR_filtered %>%
  group_by(sample_type, mechanism) %>%
  summarize(num_genes = n(), .groups = "drop") %>% 
  mutate(mechanism = factor(mechanism, levels = 
                              mechanism_ord$mechanism)) %>% 
  mutate(mechanism = fct_rev(mechanism)) %>% 
  mutate(sample_type = ifelse(sample_type == "cow", "Cow Dung", "Soil Floor")) %>% 
  mutate(mechanism = str_to_sentence(mechanism))


mech_plot <- ggplot(mechanism_by_sample_type, aes(x = sample_type, y = num_genes)) + 
  geom_bar(aes(fill = mechanism), stat = "identity", position="fill", col="black",
           width=0.8) +
  scale_fill_manual(values = qualitative_colors)+
  theme_minimal() +
  scale_y_continuous(labels = c("0","25","50","75","100"),
                     breaks = c(0,.25, .5, .75, 1)) + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black",size=12)) +
  ylab("Relative abundance of ARGs (%)")

ggsave(paste0(figure_path, "fig-AMR-mechanism_rel_abundance.pdf"),
       mech_plot, width = 7.5, height = 3)
