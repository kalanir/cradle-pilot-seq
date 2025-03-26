#!/bin/bash

R CMD BATCH 1-fig-bracken-pathogen-rel-abundance.R
R CMD BATCH 2-fig-bracken-pathogen-heatmap.R 
R CMD BATCH 3-fig-alpha-diversity-pathogen-species.R &
R CMD BATCH 4-fig-alpha-diversity-pathogen-genus.R &
R CMD BATCH 5-fig-beta-diversity-pathogen-species.R &
R CMD BATCH table-all-orgs.R

sleep 3

rm Rplots.pdf
rm .RData