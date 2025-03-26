#!/bin/bash

R CMD BATCH 1-process-AMR.R 
R CMD BATCH 2-fig-AMR-heatmap.R &
R CMD BATCH 3-fig-AMR-rel_abundance.R &
R CMD BATCH 4-fig-AMR-ARG-risk.R &
R CMD BATCH 5-fig-AMR-heatmap-risk.R

sleep 3
rm Rplots.pdf
rm .RData