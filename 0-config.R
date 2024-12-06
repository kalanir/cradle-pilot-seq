#-------------------------------------
# floor rct

# configure data directories
# source base functions
# load libraries
#-------------------------------------

library(tidyverse)
library(here)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(assertthat)
library(boxr)

# FILL THE DATA DIRECTORY
# if(Sys.getenv("USER") == ""){ 
#   box_path = "/Users/"
# }

# Define directories
proj_dir = paste0(here::here(), "/")
data_dir = here::here("data/")
figure_path = here::here("figures/")
table_path = here::here("tables/")

