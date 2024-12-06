#######################################
# CRADLE trial

# Analysis of sequencing data from 
# 2023 environmental pilot

# Descriptive summary for text  
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

species_all = readRDS(paste0(data_dir, "Bracken_processed_data/all/species_percentages.RDS"))
rownames(species_all) = species_all$name 
species_all = species_all[,-c(1:2)]

genera_all = readRDS(paste0(data_dir, "Bracken_processed_data/all/genus_percentages.RDS"))
rownames(genera_all) = genera_all$name
genera_all = genera_all[,-c(1:2)]

species_pathogen = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens/pathogen_species_percentages.RDS"))
rownames(species_pathogen) = species_pathogen$name 
species_pathogen = species_pathogen[,-c(1:2)]

genera_pathogen = readRDS(paste0(data_dir, "Bracken_processed_data/filtered_pathogens/pathogen_genus_percentages.RDS"))
rownames(genera_pathogen) = genera_pathogen$name
genera_pathogen = genera_pathogen[,-c(1:2)]



# check that percentages >0 for all species/genera
assert_that(all(rowSums(species_all))>0)
assert_that(all(rowSums(genera_all))>0)
assert_that(all(rowSums(species_pathogen))>0)
assert_that(all(rowSums(genera_pathogen))>0)


# all microbes -------------------------------------------------
## total species by sample type -------------------------------------------------
# cow
cow_species <- species_all[, grep("C", colnames(species_all))]  %>% 
  as.data.frame() %>% 
  mutate(total_species = rowSums(.)) %>% 
  filter(total_species > 0) 

nrow(cow_species)

# soil
# cow
soil_species <- species_all[, grep("S", colnames(species_all))]  %>% 
  as.data.frame() %>% 
  mutate(total_species = rowSums(.)) %>% 
  filter(total_species > 0) 

nrow(soil_species)

## total genera by sample type -------------------------------------------------
# cow
cow_genera <- genera_all[, grep("C", colnames(genera_all))]  %>% 
  as.data.frame() %>% 
  mutate(total_species = rowSums(.)) %>% 
  filter(total_species > 0) 

nrow(cow_genera)

# soil
# cow
soil_genera <- genera_all[, grep("S", colnames(genera_all))]  %>% 
  as.data.frame() %>% 
  mutate(total_species = rowSums(.)) %>% 
  filter(total_species > 0) 

nrow(soil_genera)



# pathogens -------------------------------------------------
## total species by sample type -------------------------------------------------
# cow
cow_species <- species_pathogen[, grep("C", colnames(species_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(total_species = rowSums(.)) %>% 
  filter(total_species > 0) 

nrow(cow_species)
  
# soil
# cow
soil_species <- species_pathogen[, grep("S", colnames(species_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(total_species = rowSums(.)) %>% 
  filter(total_species > 0) 

nrow(soil_species)

## total genera by sample type -------------------------------------------------
# cow
cow_genera <- genera_pathogen[, grep("C", colnames(genera_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(total_species = rowSums(.)) %>% 
  filter(total_species > 0) 

nrow(cow_genera)

# soil
# cow
soil_genera <- genera_pathogen[, grep("S", colnames(genera_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(total_species = rowSums(.)) %>% 
  filter(total_species > 0) 

nrow(soil_genera)

# average no. species per sample by sample type -------------------------------------------------
## all microbes -------------------------------------------------
colSums(species_all[, grep("C", colnames(species_all))] > 0) %>% mean()
colSums(species_all[, grep("S", colnames(species_all))] > 0) %>% mean()

colSums(genera_all[, grep("C", colnames(genera_all))] > 0) %>% mean()
colSums(genera_all[, grep("S", colnames(genera_all))] > 0) %>% mean()


## pathogens -------------------------------------------------
colSums(species_pathogen[, grep("C", colnames(species_pathogen))] > 0) %>% mean()
colSums(species_pathogen[, grep("S", colnames(species_pathogen))] > 0) %>% mean()

colSums(genera_pathogen[, grep("C", colnames(genera_pathogen))] > 0) %>% mean()
colSums(genera_pathogen[, grep("S", colnames(genera_pathogen))] > 0) %>% mean()


# most common genera  -------------------------------------------------
## all microbes -------------------------------------------------
top_cow_genera <- genera_all[, grep("C", colnames(genera_all))]  %>% 
  as.data.frame() %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance)) %>% 
  slice(1:5) 

top_soil_genera <- genera_all[, grep("S", colnames(genera_all))]  %>% 
  as.data.frame() %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance)) %>% 
  slice(1:5) 

## pathogens -------------------------------------------------
top_cow_genera_pathogens <- genera_pathogen[, grep("C", colnames(genera_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(mean_abundance = rowMeans(.),
         tot_abundance = rowSums(.)) %>% 
  arrange(desc(mean_abundance)) %>% 
  slice(1:20) 

top_soil_genera_pathogens <- genera_pathogen[, grep("S", colnames(genera_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance)) %>% 
  slice(1:10) 


# most common species  -------------------------------------------------
## pathogens -------------------------------------------------
top_cow_species_pathogens <- species_pathogen[, grep("C", colnames(species_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(mean_abundance = rowMeans(.),
         tot_abundance = rowSums(.)) %>% 
  arrange(desc(mean_abundance)) %>% 
  slice(1:10) 

top_soil_genera_pathogens <- species_pathogen[, grep("S", colnames(species_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance)) %>% 
  slice(1:10) 

# paired samples with relative abundance>10% in both samples  -------------------------------------------------
paired_common_species <- species_pathogen[,  colnames(species_pathogen) %in% c("S1","C1")]  %>% 
  as.data.frame() %>% 
  filter(S1>2 & C1>2) 

get_common_shared_species <- function(hhid){
  ids = c(paste0("S",hhid), paste0("C",hhid))
  d = species_pathogen[,  colnames(species_pathogen) %in% ids]  %>% 
    as.data.frame() 
  colnames(d) = c("S","C")
  filtered = d %>%  filter(S>5 & C>5) 
  if(nrow(filtered)>0){
    return(data.frame(hhid = hhid,
                      species = paste0(rownames(filtered))))
  }else{
    return(data.frame(hhid = hhid,
                      species = "None"))
  }

}

paired_common_species = lapply(as.list(1:10), get_common_shared_species) %>% bind_rows() 
table(paired_common_species$species)

paired_common_species %>% group_by(species) %>% summarise(nhh = length(unique(hhid)))

# pathogen species in soil but not cow dung  -------------------------------------------------
species_pathogen[, grep("S", colnames(species_pathogen))]  %>% 
  as.data.frame() %>% 
  filter(rowSums(.)>0) %>% 
  rownames() %>% 
  setdiff(
    species_pathogen[, grep("C", colnames(species_pathogen))]  %>% 
      as.data.frame() %>% 
      filter(rowSums(.)>0) %>% 
      rownames()
  ) %>% 
  length()

# pathogen species in cow dung but not soil  -------------------------------------------------
species_pathogen[, grep("C", colnames(species_pathogen))]  %>% 
  as.data.frame() %>% 
  filter(rowSums(.)>0) %>% 
  rownames() %>% 
  setdiff(
    species_pathogen[, grep("S", colnames(species_pathogen))]  %>% 
      as.data.frame() %>% 
      filter(rowSums(.)>0) %>% 
      rownames()
  ) %>% 
  length()

# pathogen species in both cow dung and soil  -------------------------------------------------
species_pathogen[, grep("C", colnames(species_pathogen))]  %>% 
  as.data.frame() %>% 
  filter(rowSums(.)>0) %>% 
  rownames() %>% 
  intersect(
    species_pathogen[, grep("S", colnames(species_pathogen))]  %>% 
      as.data.frame() %>% 
      filter(rowSums(.)>0) %>% 
      rownames()
  ) %>% 
  length()

# Pathogens that were present in at least 5 of 10 households in both soil and cow dung -------------------------------------------------
half_hh_species_cow <- species_pathogen[, grep("C", colnames(species_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(mean_abund =apply(., 1, mean)) %>% 
  filter(mean_abund>1) %>% 
  dplyr::select(-mean_abund) %>% 
  mutate(at_least_5_nonzero = rowSums(. > 0) >= 5) %>% 
  filter(at_least_5_nonzero) %>% 
  rownames()
  
half_hh_species_soil <- species_pathogen[, grep("S", colnames(species_pathogen))]  %>% 
  as.data.frame() %>% 
  mutate(mean_abund =apply(., 1, mean)) %>% 
  filter(mean_abund>1) %>% 
  dplyr::select(-mean_abund) %>% 
  mutate(at_least_5_nonzero = rowSums(. > 0) >= 5) %>% 
  filter(at_least_5_nonzero) %>% 
  rownames()

intersect(half_hh_species_cow, half_hh_species_soil) 
