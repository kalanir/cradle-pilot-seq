#######################################
# CRADLE trial

# Visualization of processed sequencing AMR 
# data from 2023 environmental pilot
#######################################
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))


AMR_reports_ids = readRDS(paste0(data_dir,
                                 "CZ_processed_data/AMR/AMR_data_raw.RDS"))
AMR_filtered <- AMR_reports_ids %>% filter((reads >=5 & read_coverage_breadth >= 10) |
                                             contig_coverage_breadth >=10)

# ----------------------------------------------------------------
# Langelier: https://pmc.ncbi.nlm.nih.gov/articles/PMC11092797/
# If in the same sample one AMR gene was found by the read approach and a different 
# AMR gene from the same gene family was found by the contig approach, the first AMR 
# gene was omitted and only the second AMR gene was plotted. 
# Rationale: sometimes short reads alone cannot sufficiently distinguish between 
# highly similar alleles or genes from the same gene family.
# ----------------------------------------------------------------
# Create a new variable 'keep' initialized with TRUE for all rows
data <- AMR_filtered
data$keep <- TRUE

# Get the unique hhid values
unique_sampleids <- unique(data$sample_id)

# Iterate over each unique hhid
for (sampleid in unique_sampleids) {
  # Subset the data for the current hhid
  sampleid_data <- data[data$sample_id == sampleid, ]
  
  # Get the unique gene_family values within the current hhid
  unique_gene_families <- unique(sampleid_data$gene_family)
  
  # Iterate over each unique gene_family within the current hhid
  for (gene_family in unique_gene_families) {
    # Subset the data for the current gene_family within the current hhid
    gene_family_data <- sampleid_data[sampleid_data$gene_family == gene_family, ]
    
    # Check if there are multiple genes in the same gene_family
    if (nrow(gene_family_data) > 1) {
      # Check if there's a gene found by the contig approach (contigs > 0)
      contig_gene <- gene_family_data[gene_family_data$contigs > 0, ]
      
      if (nrow(contig_gene) > 0) {
        # If a contig gene is found, set 'keep' to FALSE for all read-based genes (contigs == 0)
        gene_family_data$keep[gene_family_data$contigs == 0] <- FALSE
      }
      
    }
    
    # Update the 'keep' variable in the original 'data' data frame
    data[data$sample_id == sampleid & data$gene_family == gene_family, "keep"] <- gene_family_data$keep
    
  }
}

# Subset the data to keep only the rows where 'keep' is TRUE
filtered_data <- data[data$keep, ]

saveRDS(filtered_data, paste0(data_dir,
                              "/CZ_processed_data/AMR/AMR_data_filtered.RDS"))
