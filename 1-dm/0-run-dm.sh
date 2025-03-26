#!/bin/bash

R CMD BATCH 1-extract_czi_pathogen_list.R
sh 2-download-rawdata.sh