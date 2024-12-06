#######################################
# CRADLE trial

# Visualization of processed sequencing AMR 
# data from 2023 environmental pilot

# Plot ARG risk heatmap
#######################################
rm(list=ls())
library(readxl)
source(paste0(here::here(), "/0-config.R"))

## download supp file from Zhang et al., 2022 https://doi-org.stanford.idm.oclc.org/10.1038/s41467-022-29283-8
download.file('https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-29283-8/MediaObjects/41467_2022_29283_MOESM13_ESM.xlsx', destfile = paste0(data_dir,"Zhang 2022 supplement data 10.xlsx"), method = "wget")
zhang = read_excel(paste0(data_dir,"Zhang 2022 supplement data 10.xlsx"))


zhang = zhang %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "Nocardia rifampin resistant beta-subunit of RNA polymerase (rpoB2)", "rpoB2", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "tetQ", "tet(Q)", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "tetW", "tet(W)", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "Bifidobacterium adolescentis rpoB conferring resistance to rifampicin", "rpoB", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "AAC(6')-Ie-APH(2'')-Ia", "AAC(6')-Ie-APH(2'')-Ia bifunctional protein", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "tetO", "tet(O)", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "tet32", "tet(32)", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "Bifidobacteria intrinsic ileS conferring resistance to mupirocin", "ileS", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "mexK", "MexK", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "tet44", "tet(44)", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "mexW", "MexW", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "Tet(X4)", "tet(X4)", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "acrF", "AcrF", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "yojI", "YojI", `ARO Term`)) %>%
  mutate(`ARO Term` = ifelse(`ARO Term` == "mexI", "MexI", `ARO Term`)) 

AMR_gene_data = readRDS(paste0(data_dir, 
                              "/CZ_processed_data/AMR/AMR_data_filtered.RDS"))
  
merged = left_join(AMR_gene_data, zhang, by = c("gene" ="ARO Term"))
table(is.na(merged$Rank))

table(merged$Rank, merged$sample_id)

df= as.data.frame(table(merged$Rank, merged$sample_id)) %>% 
  rename(risk = Var1,
         sample_id = Var2,
         ARGcount = Freq) %>% 
  mutate(sample_type = ifelse(sample_id %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"), 
  "Cow dung", "Soil floor")) %>% 
  mutate(risk = as.character(risk)) %>%
  mutate(risk = ifelse(risk == "RI=0", "No risk", risk)) %>% 
  mutate(risk = factor(risk, levels = c("No risk", "Q4", "Q3", "Q2", "Q1"))) %>% 
  mutate(sample_id = factor(sample_id, levels = c(
    "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"
  )
  ))
  
# make heat map
plot <- ggplot(df, aes(x = sample_id, y = risk, fill = ARGcount)) +
  geom_tile() +
  scale_fill_viridis(option="B", end = 0.8) +
  theme_minimal() +
  facet_grid(~sample_type, scales = "free") + 
  labs(x = "Sample ID",
       y = "Risk category",
       fill = "ARG count") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=11))

plot

ggsave(plot, filename = paste0(figure_path, "fig-AMR-risk.pdf"),width = 8, height = 2.25, dpi = 300)
