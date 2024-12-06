#######################################
# CRADLE trial

# Visualization of processed sequencing AMR 
# data from 2023 environmental pilot
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# ----------------------------------------------------------------
# ID correction
# ----------------------------------------------------------------
AMR_reports = readRDS(paste0(data_dir,"/CZ_processed_data/AMR/AMR_data_raw.RDS"))
corrected_ids <- readRDS(paste0(data_dir,"/CZ_processed_data/sequenced_ids.RDS"))

AMR_reports_ids <- left_join(AMR_reports, corrected_ids %>% dplyr::select(-c(Sample_id_F, Sample_id_MSC)),
                             by = "sample_name") %>% 
  mutate(sample_id = ifelse(grepl("C", sample_name), paste0("C", hhid), paste0("S", hhid))) %>% 
  dplyr::select(-sample_name) 


AMR_filtered <- AMR_reports_ids %>% filter(reads >=5 & (read_coverage_breadth >= 10 |
                                                      contig_coverage_breadth >=10))

# ----------------------------------------------------------------
# Only keep genes with contigs == 0 if there is only one gene in the gene family
# ----------------------------------------------------------------
# Create a new variable 'keep' initialized with TRUE for all rows
data <- AMR_filtered
data$keep <- TRUE

# Get the unique hhid values
unique_hhids <- unique(data$hhid)

# Iterate over each unique hhid
for (hhid in unique_hhids) {
  # Subset the data for the current hhid
  hhid_data <- data[data$hhid == hhid, ]
  
  # Get the unique gene_family values within the current hhid
  unique_gene_families <- unique(hhid_data$gene_family)
  
  # Iterate over each unique gene_family within the current hhid
  for (gene_family in unique_gene_families) {
    # Subset the data for the current gene_family within the current hhid
    gene_family_data <- hhid_data[hhid_data$gene_family == gene_family, ]
    
    # Check if there are multiple genes in the same gene_family
    if (nrow(gene_family_data) > 1) {
      # Count the number of genes with contig == 0
      count_contig_zero <- sum(gene_family_data$contigs == 0)
      
      # If there is only one gene with contig == 0, set 'keep' to FALSE for that gene
      if (count_contig_zero == 1) {
        gene_family_data$keep[gene_family_data$contigs == 0] <- FALSE
      }
      # If there are multiple genes with contig == 0, set 'keep' to FALSE for genes with contig > 0
      else if (count_contig_zero > 1) {
        gene_family_data$keep[gene_family_data$contigs > 0] <- FALSE
      }
    }
    
    # Update the 'keep' variable in the original 'data' data frame
    data[data$hhid == hhid & data$gene_family == gene_family, "keep"] <- gene_family_data$keep
  }
}

# Subset the data to keep only the rows where 'keep' is TRUE
filtered_data <- data[data$keep, ]

saveRDS(filtered_data, paste0(data_dir,"/CZ_processed_data/AMR/AMR_data_filtered.RDS"))

