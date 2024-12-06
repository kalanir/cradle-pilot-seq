#!/bin/bash
#SBATCH --job-name=metagenomics-kraken
#SBATCH --begin=now
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=100G
#SBATCH --time=2-00:00:00
#SBATCH -o ./Report/output-classification.out # STDOUT

ml python/3.9.0
ml biology

cd $GROUP_HOME/cradle-pilot-seq/

export PYTHONPATH= <path to python>

# Snakefile cloned from https://github.com/bhattlab/kraken2_classification/ into gihub root directory
snakemake \
-s $GROUP_HOME/public/kraken2_classification/Snakefile \
--configfile 2-mNGS/1-classification-pipeline/snakemake_config_files/config.yaml \
--jobs 1 





