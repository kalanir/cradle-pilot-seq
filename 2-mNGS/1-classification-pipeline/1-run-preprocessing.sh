#!/bin/bash
#SBATCH --job-name=metagenomics-preprocessing
#SBATCH --begin=now
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=50G
#SBATCH --ntasks=50
#SBATCH --time=2-00:00:00
#SBATCH -o ./Report/output-preprocessing.out # STDOUT

EXPORT HOME_PATH=<path_to_git_directory>
ml python/3.9.0
ml biology
ml bwa

cd $GROUP_HOME/cradle-pilot-seq/

export PYTHONPATH= <path to python>

# Snakefile cloned from https://github.com/bhattlab/bhattlab_workflows/
snakemake --snakefile $GROUP_HOME/public/bhattlab_workflows/preprocessing/preprocessing.snakefile \
--configfile 2-mNGS/1-classification-pipeline/snakemake_config_files/config_preprocessing.yaml \
--jobs 100 \
--use-singularity 




