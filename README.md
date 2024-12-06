# Microbiomes and Resistomes in Household Environments with Soil Floors and Domestic Animal Cohabitation: A Study in Rural Bangladesh

## Overview
In low-income countries, inadequate housing (e.g., soil floors) and animal cohabitation may expose children to fecal organisms, increasing their risk of enteric and antimicrobial-resistant infections. Our objective was to understand whether cow cohabitation in homes with soil floors in rural Bangladesh contributed to pathogen and antimicrobial resistance genes (ARGs) in the home. Here, across 10 randomly selected households in rural Sirajganj District, we performed shotgun metagenomic sequencing of sampled floor soil and cow dung, and analyzed resulting sequencing data for pathogens and ARGs. 


## Additional Information
To run the scripts, the data directory for the user must be changed in **`0-config.R`**. This will allow for replication of study findings using scripts in this repository. Similarly, directory statement changes will be needed wherever output files are saved down. In particular, raw fastq data will need to be downloaded from NCBI Bioproject: PRJNA1130536, thereby paths must be inputed in `2-download-rawdata.sh`. We will update this once the project becomes public.

## Directory structure

**`0-analysis-bash.sh` :** bash script, replicates the entire projectfdsa`

**`0-config.R` :** configuration file that sets data directories, sources base functions, and loads required libraries

**`1-dm` :** folder containing data management scripts. To rerun all scripts in this subdirectory, run the bash script `0-run-dm.sh`.    
* `1-extract_czi_pathogen_list.R` : generate pathogen csv from the HTML output from https://czid.org/pathogen_list/
* `2-download-rawdata.sh` : download raw fastq data into raw_data/ folder. Requires paths to be inputed `<path to fastq>` 
* `3-download-reference-genome.sh`: download cow reference genome

**`2-mNGS` :** folder containing processing and analysis scripts for pathogen analysis. 
* `1-classification-pipeline` : Folder that includes scripts to process fastq files via kraken
    * `1-run-preprocessing.sh` : bash script to run pipeline for sequence preprocessing (see note at bottom)
    * `2-run-kraken.sh` : bash script to run pipeline for kraken/bracken classification (see note at bottom)
    * `snakemake_config_files` : Folder that includes snakemake configuration files to change and run kraken 
* `2-analysis` : Folder that contains scripts to run analyses on the kraken/bracken output
    * `1-post_classification_workflow.R` : estimates relative abundance for all detected genera/species
    * `2-filtered_pathogen_classification_workflow.R` : estimates relative abundance for pathogens
    * `3-summarize_read_counts.R ` : enumerataes the number of reads retained at each stage of pre-processing
    * `4-pathogen-descriptive.R` : get descriptions of detected pathogens for manuscript text
* `3-figure-scripts` : test unmeasured mediator-outcome confounding. To rerun all scripts in this subdirectory, run the bash script `0-fig-bash.sh`.           
    * `1-fig-bracken-pathogen-rel-abundance.R` : script to generate relative abundance plots.
    * `2-fig-bracken-pathogen-heatmap.R` : script to generate pathogen heatmap figures.
    * `3-fig-all-alpha-diversity-indexes.R` : script to generate alpha diversity figures across all microbes detected in kraken/bracken results.
    * `3-fig-alpha-diversity-indexes.R` : script to generate alpha diversity figures across pathogen-filtered kraken/bracken results.
    * `4-fig-all-beta-diversity-indexes.R` : script to generate beta diversity figures across all microbes detected in kraken/bracken results.
    * `4-fig-beta-diversity-indexes.R` : script to generate beta diversity figures across pathogen-filtered kraken/bracken results.    
    * `table-all-orgs.R` : script to generate pathogen count tables.  

**`3-AMR` :** folder containing analysis scripts for AMR analysis. To rerun all scripts in this subdirectory, run the bash script `0-run-AMR.sh`.     
* `1-process-AMR.R` : script to process and filtering raw CZID AMR output data. 
* `2-fig-AMR-heatmap.R` : script to generate AMR heatmap.   
* `3-fig-AMR-barplot.R` : script to generate AMR barplots.
* `4-fig-AMR-rel_abundance.R` : script to generate AMR relative abundance barplots.
* `5-fig-AMR-ARG-risk.R` : script to generate AMR risk heatmap.

**`data` :** folder containing table scripts. Before running scripts in this subdirectory, first run `2-analysis/0-analysis-bash.sh` to recreate the result data sets in `7-results`.    
* `pathogen_list.txt` : html import of  https://czid.org/pathogen_list/
* `czi_pathogen_list.csv` : readable pathogen csv of `pathogen_list.txt` after running `1-dm/2-download-rawdata.sh`
* `Bracken_processed_data` : folder that contains data generated by CZID
    * `id_master_list.xlsx` : xlsx file to convert sample names to associated sampleid (associated with hhid)
    * `all` : folder that contains full kraken/bracken output
        * `genus_counts.RDS` : file of genus counts from kraken/bracken output
        * `genus_percentages.RDS` : file of genus percentages from kraken/bracken output
        * `species_counts.RDS` : file of species counts from kraken/bracken output
        * `species_percentages.RDS` : file of species percentages from kraken/bracken output
    * `filtered_pathogens` : folder that contains pathogen-filtered karaken/bracken output
        * `pathogen_genus_counts.RDS` : file of pathogen-filtered genus counts from kraken/bracken output
        * `pathogen_genus_percentages.RDS` : file of pathogen-filtered genus percentages from kraken/bracken output
        * `pathogen_species_counts.RDS` : file of pathogen-filtered species counts from kraken/bracken output
        * `pathogen_species_percentages.RDS` : file of pathogen-filtered species percentages from kraken/bracken output
* `CZ_processed_data` : folder that contains data generated by CZID
    * `sequenced_ids.RDS` : rds file to convert sample names to associated sampleid (associated with hhid)
    * `AMR` : folder that contains combined data from CZID AMR pipeline
        * `AMR_data_raw.RDS` : file of the combined raw CZID AMR data across samples - by sample_name
        * `AMR_data_filtered.RDS` : file of the combined filtered CZID AMR data across samples - by sample_id

**`figures` :** folder containing figure files.

**`6-tables` :** folder containing table files.
* `table_A1.csv` : sample metadata table A1 from paper
* `TableA2_sample_readcounts.csv` : sample fastq read counts from kraken/bracken pipeline
* `table_cow_genus.csv` : output from `3-figure-scripts/table-all-orgs.R`
* `table_cow_species.csv` : output from `3-figure-scripts/table-all-orgs.R`
* `table_soil_genus.csv` : output from `3-figure-scripts/table-all-orgs.R`
* `table_soil_species.csv` : output from `3-figure-scripts/table-all-orgs.R`
* `sample_readcounts.tsv` : readcount details

Contributors: Jade Benjamin-Chung, Anna T. Nguyen, Gabby B. Heitmann, Kalani Ratnasiri

This analysis leveraged public workflows that were originally authored by members of the Bhatt Lab (Stanford Department of Genetics):
* Preprocessing: https://github.com/bhattlab/bhattlab_workflows/
* Kraken2 + Bracken Classification, Post-Classification Processing: https://github.com/bhattlab/kraken2_classification
