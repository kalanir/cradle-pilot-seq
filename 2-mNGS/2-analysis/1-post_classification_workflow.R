# Processing kraken results into matrices and plots
# Ben Siranosian - Bhatt lab - Stanford Genetics
# bsiranosian@gmail.com
# January 2019 - June 2023
proj_dir = paste0(here::here(),"/")
suppressMessages(library(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(rafalib, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(vegan, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(reshape2, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(RColorBrewer, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(cmapR, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(compositions, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(zCompositions, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(ALDEx2, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(ggpubr, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(dplyr, quietly = TRUE, warn.conflicts = FALSE))
select = dplyr::select

options(stringsAsFactors = F)

# options we need from snakemake
sample_reads_file <- paste0(proj_dir,'2-mNGS/1-classification-pipeline/snakemake_config_files/classification_input.txt')
sample_groups_file <- paste0(proj_dir,'2-mNGS/1-classification-pipeline/snakemake_config_files/sample_groups.tsv')
workflow.outdir <-  proj_dir

result.dir <-  paste0(proj_dir,'results/')
use.bracken.report <- TRUE
classification_method = ifelse(use.bracken.report, "bracken", "kraken")
if (classification_method == "kraken") {
  result.dir = paste0(result.dir, "-kraken")
}
remove.chordata <- TRUE

# Use locations in the working directory for scripts and taxonomy array
scripts.folder <- file.path(workflow.outdir, "scripts")
tax.array.file <- file.path(workflow.outdir, "taxonomy_array.tsv")

# Downstream processing for filtering OTUs
## This removes any classification result where all samples are below this threshold,
## which serves to remove the very long tail of lowly abundant species in Kraken2/Bracken results.
## This parameter can be tuned, but I recommend you keep something here to reduce the long tail.
min_otu_percentage <- 0

# Set up directories and create those that don't exist
classification.folder <- file.path(workflow.outdir, 'classification')

# Load other data processing and plotting scripts
source.script.process <- file.path(scripts.folder, 'process_classification_gctx.R')
suppressMessages(source(source.script.process))

# read from sample_reads_file or sample_reports_file as input
if(sample_reads_file != ""){
  sample.reads <- read.table(sample_reads_file, sep='\t', quote='', header=F, comment.char = "#", colClasses = 'character')
  colnames(sample.reads) <- c('sample', 'r1', 'r2')[1:ncol(sample.reads)]
  # ensure we skip first row if it's a comment
  if ((tolower(sample.reads[1,1]) == 'sample') | (substr(sample.reads[1,1], 1,1) == "#")){
    sample.reads <- sample.reads[2:nrow(sample.reads), ]
  }
  # create dataframe of output files that would be in sample.reports
  sample.reports <- data.frame(sample = sample.reads$sample, 
                               kraken_report = sapply(sample.reads$sample, function(x) file.path(classification.folder, paste(x, '.krak.report', sep=''))),
                               bracken_report = sapply(sample.reads$sample, function(x) file.path(classification.folder, paste(x, '.krak_bracken_species.report', sep=''))))                    
} else if(sample_reports_file != ""){
  sample.reports <- read.table(sample_reports_file, sep='\t', quote='', header=F, comment.char = "#", colClasses = 'character')
  colnames(sample.reports) <- c('sample', 'kraken_report', 'bracken_report')[1:ncol(sample.reports)]
  # ensure we skip first row if it's a comment
  if ((tolower(sample.reports[1,1]) == 'sample') | (substr(sample.reports[1,1], 1,1) == "#")){
    sample.reports <- sample.reports[2:nrow(sample.reports), ]
  }
}
sample.names <- sample.reports$sample
sample.number <- nrow(sample.reports)

# Define files for loading in based on sample.
if(use.bracken.report){
  flist <- sample.reports[,3]
} else {
  flist <- sample.reports[,2]
}
names(flist) <- sample.names

# Ensure all expected results files exist
if (!(all(file.exists(flist)))){
  # print which files don't exist
  message('The follwing files do not exist:')
  print(flist[!sapply(flist, file.exists)])
  stop("Some classification files do not exist!")
}

# if a groups file is specified, read it. Otherwise assign everything to one group
if (sample_groups_file != '') {
  sample.groups <- read.table(sample_groups_file, sep='\t', quote='', header=F, comment.char = "#", colClasses = 'character')
  colnames(sample.groups) <- c('sample', 'group')
  # if first sample is sample, discard the line
  if(tolower(sample.groups[1,'sample']) == 'sample'){
    sample.groups <- sample.groups[2:nrow(sample.groups),]
  }
} else {
  sample.groups <- data.frame(sample=sample.names, group='All')
}

# get taxonomy array
# classification at each level for each entry in the database.
# simplified after processing krakens default taxonomy
tax.array <- read.table(tax.array.file, sep='\t', quote='', header=F, comment.char = '', colClasses = 'character')
colnames(tax.array) <- c('id', 'taxid', 'root', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies')[1:ncol(tax.array)]
# Bug in generation code gave muliple zero taxids. Eliminate that here now
tax.array <- tax.array[!duplicated(tax.array$taxid),]
rownames(tax.array) <- tax.array$taxid
# some duplicate names in this. if so, change them to include the tax level
dup.ids <- tax.array$id[duplicated(tax.array$id)]
dup.inds <- which(tax.array$id %in% dup.ids)
# fix these by adding taxid to the name
tax.array[tax.array$id %in% dup.ids, "id"] <-
  paste(tax.array[tax.array$id %in% dup.ids, "id"], ' (', tax.array[tax.array$id %in% dup.ids, "taxid"], ')', sep='')

# ensure reads and groups have the same data
if (!(all(sample.groups$sample %in% sample.names) & all(sample.names %in% sample.groups$sample))){
  message('sample.reads:')
  print(sample.reads)
  message('sample.groups:')
  print(sample.groups)
  message('sample.groups$sample[!(sample.groups$sample %in% sample.names)]')
  print(sample.groups$sample[!(sample.groups$sample %in% sample.names)])
  message('sample.names[!(sample.names %in% sample.groups$sample)]')
  print(sample.names[!(sample.names %in% sample.groups$sample)])
  stop('Sample reads and sample groups dont contain the same samples... check inputs')
}

# special case for UHGG and MAG databases.
## If reading from UHGG, the taxonomy levels go R, R1-R7
test.df <- kraken_file_to_df(flist[1])
## UHGG database will have this structure
uhgg <- all(paste0('R', 1:7) %in% test.df$tax.level)
# MAG database will have this structure
## Number of genus classifications is zero
segata <- !uhgg & sum(sum(test.df$tax.level=='G')) == 0

# load classification results from each sample and process into gct format
message(paste('Loading data from ', length(flist), ' kraken/bracken results.', sep=''))
df.list <- lapply(flist, function(x) kraken_file_to_df(x))
# merge classification matrix
merge.mat <- merge_kraken_df_list(df.list)
# sample metadata just has groups for now
sample.metadata <- data.frame(id=sample.groups$sample, group=sample.groups$group)
# construct gct object from reads matrix, sample metadata, row metadata
kgct <- make_gct_from_kraken(merge.mat, sample.metadata, tax.array)

# remove Chordata reads if desired
if(remove.chordata){
  kgct <- subset_gct(kgct, rid=kgct@rid[kgct@rdesc$phylum != "Chordata"])
}
# filter to each of the taxonomy levels
if (segata){
  filter.levels <- c('species')
  kgct.filtered.list <- list(species=subset_gct(kgct, rid=kgct@rid[kgct@rid != 'root']))
} else {
  filter.levels <-  c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
  kgct.filtered.list <- lapply(filter.levels, function(level) {
    subset_kgct_to_level(kgct, level)
  })
  names(kgct.filtered.list) <- filter.levels
}
# print("kgct.filtered.list[['species']]@mat")
# print(kgct.filtered.list[['species']]@mat)

# subset to classified taxa only
unclassified.rownames <- c('unclassified', 'classified at a higher level')
kgct.filtered.classified.list <- lapply(kgct.filtered.list, function(x) subset_gct(x, rid=x@rid[!(x@rid %in% unclassified.rownames)]))

# normalize to percentages
kgct.filtered.percentage.list <- lapply(kgct.filtered.list, function(x) normalize_kgct(x, min_otu_percentage=min_otu_percentage))
kgct.filtered.classified.percentage.list <- lapply(kgct.filtered.classified.list, function(x) normalize_kgct(x, min_otu_percentage=min_otu_percentage))

extract_taxid = function(classified_list){
  classified_list@mat %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("name") %>% 
    bind_cols(classified_list@rdesc %>% tibble::rownames_to_column("name") %>% select(taxid)) %>% 
    select(name, taxid, everything())
}


saveRDS(kgct.filtered.classified.list$species %>% extract_taxid(), 
        here::here("data","Bracken_processed_data", "contigs","/species_counts.RDS"))
saveRDS(kgct.filtered.classified.list$genus %>% extract_taxid(), 
        here::here("data","Bracken_processed_data", "contigs","/genus_counts.RDS"))

saveRDS(kgct.filtered.classified.percentage.list$species %>% extract_taxid(), 
        here::here("data","Bracken_processed_data", "contigs","/species_classified_percentages.RDS"))
saveRDS(kgct.filtered.classified.percentage.list$genus %>% extract_taxid(), 
        here::here("data","Bracken_processed_data", "contigs","/genus_classified_percentages.RDS"))

