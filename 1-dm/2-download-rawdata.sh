#!/bin/bash
EXPORT HOME_PATH=<path_to_git_directory>
mkdir -p $HOME_PATH/cradle-metagenomic-classification/raw_data/
cd $HOME_PATH/cradle-metagenomic-classification/

curl -L <path to fastq> --output raw_data/C10_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C10_R2.fastq.gz
curl -L <path to fastq> --output raw_data/C1_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C1_R2.fastq.gz
curl -L <path to fastq> --output raw_data/C2_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C2_R2.fastq.gz
curl -L <path to fastq> --output raw_data/C3_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C3_R2.fastq.gz
curl -L <path to fastq> --output raw_data/C4_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C4_R2.fastq.gz
curl -L <path to fastq> --output raw_data/C5_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C5_R2.fastq.gz
curl -L <path to fastq> --output raw_data/C6_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C6_R2.fastq.gz
curl -L <path to fastq> --output raw_data/C7_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C7_R2.fastq.gz
curl -L <path to fastq> --output raw_data/C8_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C8_R2.fastq.gz
curl -L <path to fastq> --output raw_data/C9_R1.fastq.gz
curl -L <path to fastq> --output raw_data/C9_R2.fastq.gz
curl -L <path to fastq> --output raw_data/S10_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S10_R2.fastq.gz
curl -L <path to fastq> --output raw_data/S1_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S1_R2.fastq.gz
curl -L <path to fastq> --output raw_data/S2_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S2_R2.fastq.gz
curl -L <path to fastq> --output raw_data/S3_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S3_R2.fastq.gz
curl -L <path to fastq> --output raw_data/S4_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S4_R2.fastq.gz
curl -L <path to fastq> --output raw_data/S5_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S5_R2.fastq.gz
curl -L <path to fastq> --output raw_data/S6_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S6_R2.fastq.gz
curl -L <path to fastq> --output raw_data/S7_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S7_R2.fastq.gz
curl -L <path to fastq>  --output raw_data/S8_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S8_R2.fastq.gz
curl -L <path to fastq> --output raw_data/S9_R1.fastq.gz
curl -L <path to fastq> --output raw_data/S9_R2.fastq.gz
  