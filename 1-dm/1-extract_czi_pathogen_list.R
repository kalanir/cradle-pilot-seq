library(stringr)
library(unglue)
library(tidyverse)

# Read in HTML output from https://czid.org/pathogen_list/
pathogen_list = read.delim(here::here("data", "pathogen_list.txt"))
pathogen_list = paste0(pathogen_list$X.html..head., collapse = "\n")

split1 = str_split(pathogen_list, "<h2>Bacteria</h2>")[[1]]
split2 = str_split(split1[2], "<h2>Eukaryota</h2>")[[1]]
split3 = str_split(split2[2], "<h2>Viruses</h2>")[[1]]

bacteria_list = split2[1]
eukaryota_list = split3[1]
virus_list = split3[2]

extract_path_df = function(path_list, path_type) {
  split_pathogens = str_split(path_list, 'pathogenName')[[1]]
  split_pathogens = split_pathogens[2:length(split_pathogens)]
  path_df = unglue_data(split_pathogens, "-1eRM7>{name}</div><a href={ncbi_link} {=.*}Tax ID: {tax_id}</a>{=.*}") %>% 
    mutate(type = path_type) %>% 
    select(type, name, tax_id, ncbi_link)
  return (path_df)
}

path_df = mapply(extract_path_df, 
                 path_list = c(bacteria_list, eukaryota_list, virus_list), 
                 path_type = c("Bacteria", "Eukaryote", "Virus"), 
                 SIMPLIFY = F) %>% 
  bind_rows()

write_csv(path_df, here::here("data", "czi_pathogen_list.csv"))
