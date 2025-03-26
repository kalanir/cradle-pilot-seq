#!/bin/bash
EXPORT HOME_PATH=<path_to_git_directory>
mkdir -p $HOME_PATH/cradle-metagenomic-classification/raw_data/
cd $HOME_PATH/cradle-metagenomic-classification/

curl -L <path to fastq> --output raw_data/C10_R1.fastq.gz
curl -L <path to fastq> --output data/contigs/C1_2_533733_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/C2_1_533734_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/C3_1_533735_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/C4_1_535072_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/C5_1_535073_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/C6_1_549008_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/C7_1_549010_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/C8_1_549011_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/C9_1_549017_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/C10_1_549021_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S1_506125_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S2_506124_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S3_506123_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S4_1_506251_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S5_1_506252_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S6_1_506253_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S7_506752_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S8_507214_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S9_507302_contigs_nh.fasta
curl -L <path to fastq> --output data/contigs/S10_507395_contigs_nh.fasta
  