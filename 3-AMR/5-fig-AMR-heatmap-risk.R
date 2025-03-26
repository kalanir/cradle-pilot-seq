#######################################
# CRADLE trial

# Visualization of processed sequencing AMR 
# data from 2023 environmental pilot

# Plot ARG risk heatmap
#######################################
rm(list=ls())
library(readxl)
library(ggh4x)
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
  mutate(`ARO Term` = ifelse(`ARO Term` == "mexI", "MexI", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanWB", "vanW gene in vanB cluster", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanWI", "vanW gene in vanI cluster", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanYA", "vanY gene in vanA cluster", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanYB", "vanY gene in vanB cluster", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanYF", "vanY gene in vanF cluster", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanYG", "vanY gene in vanG cluster", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanYM", "vanY gene in vanM cluster", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanWG", "vanW gene in vanG cluster", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanTG", "vanT gene in vanG cluster", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "tolC", "TolC", `ARO Term`)) %>%  
  mutate(`ARO Term` = ifelse(`ARO Term` == "vanRO", "vanR gene in vanO cluster", `ARO Term`)) 

AMR_gene_data = readRDS(paste0(data_dir, 
                               "CZ_processed_data/AMR/AMR_data_filtered.RDS")) 

merged = left_join(AMR_gene_data, zhang, by = c("gene" ="ARO Term"))
table(is.na(merged$Rank))

table(merged$Rank, merged$sample_id)


# gene heatmap ------------------------------------------------------------------
plot_data = merged %>% 
  mutate(gene = as.character(gene)) %>% 
  mutate(gene = ifelse(gene=="vanY gene in vanM cluster","vanYM",gene)) %>% 
  mutate(gene = ifelse(gene=="vanY gene in vanG cluster","vanYG",gene)) %>% 
  mutate(gene = ifelse(gene=="vanY gene in vanF cluster","vanYF",gene)) %>% 
  mutate(gene = ifelse(gene=="vanY gene in vanB cluster","vanYB",gene)) %>% 
  mutate(gene = ifelse(gene=="vanY gene in vanA cluster","vanYA",gene)) %>% 
  mutate(gene = ifelse(gene=="vanXY gene in vanG cluster","vanXYG",gene)) %>% 
  mutate(gene = ifelse(gene=="vanW gene in vanI cluster","vanWI",gene)) %>% 
  mutate(gene = ifelse(gene=="vanW gene in vanG cluster","vanWG",gene)) %>% 
  mutate(gene = ifelse(gene=="vanW gene in vanB cluster","vanWB",gene)) %>% 
  mutate(gene = ifelse(gene=="vanT gene in vanG cluster","vanTG",gene)) %>% 
  mutate(gene = ifelse(gene=="Neisseria gonorrhoeae 23S rRNA with mutation conferring resistance to azithromycin","Ngon23SAZM",gene)) %>% 
  mutate(gene = ifelse(gene=="Mycobacteroides abscessus 23S rRNA with mutation conferring resistance to clarithromycin","Mabs23SCLR",gene)) %>% 
  mutate(gene = ifelse(gene=="Escherichia coli soxR with mutation conferring antibiotic resistance","EcolsoxRMULT",gene)) %>% 
  mutate(gene = ifelse(gene=="Clostridioides difficile 23S rRNA with mutation conferring resistance to erythromycin and clindamycin","Cdif23SMULT",gene)) %>% 
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
  mutate(gene = as.factor(gene))   %>% 
  mutate(Rank = as.character(Rank)) %>% 
  mutate(Rank = case_when(
    Rank == "RI=0" ~ "Risk index=0",
    Rank == "Q4" ~ "Quartile 4",
    Rank == "Q3" ~ "Quartile 3",
    Rank == "Q2" ~ "Quartile 2",
    Rank == "Q1" ~ "Quartile 1")) %>% 
  filter(!is.na(Rank))



pdf(paste0(figure_path, "fig-AMR-heatmap-risk.pdf"), width=7, height=10.5)
ggplot(plot_data, aes(x = sample_id, y = gene)) + 
  geom_tile(aes(fill=Rank), col="black") + 
  scale_fill_manual("ARG Risk", values = c("#e63946","#ff794c","#ffbf5e","#ffeeba","grey")) +
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(color="white"),
    strip.background = element_rect(color="black")
  )
dev.off()

