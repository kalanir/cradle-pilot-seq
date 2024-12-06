#!/bin/bash

R CMD BATCH 1-fig-bracken-pathogen-rel-abundance.R
R CMD BATCH 2-fig-bracken-pathogen-heatmap.R 
R CMD BATCH 3-fig-all-alpha-diversity-indexes.R &
R CMD BATCH 3-fig-alpha-diversity-indexes.R &
R CMD BATCH 4-fig-all-beta-diversity-indexes.R &
R CMD BATCH 4-fig-beta-diversity-indexes.R 
R CMD BATCH table-all-orgs.R

sleep 3

rm Rplots.pdf
rm .RData