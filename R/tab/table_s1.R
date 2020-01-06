#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/tab/table_s1.R metadata/sample_info.txt metadata/new_sample_accession_nrs.txt
# ===============================================================
# This script produces Table X1 of the study "The genomic origins of a marine radiation"
# by Hench, McMillan an Puebla
#   ---------------------------------------------------------------
# ===============================================================
# args <- c('metadata/sample_info.txt', 'metadata/new_sample_accession_nrs.txt')
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)

# config -----------------------
sample_metadata_file <- as.character(args[1])
accession_file <- as.character(args[2])
# start script -------------------

sample_metadata <- read_tsv(sample_metadata_file)
accession_data <- read_tsv(accession_file)

table_out <- accession_data %>% 
  left_join(sample_metadata) %>%
  rename(ID = 'id',
         Label = 'label',
         Species = 'spec',
         Location = 'geo',
         Date = 'date',
         Latitude = 'coord_N',
         Longitude = 'coord_W',
         `Accession Number` = 'ena_accession') %>%
  select(ID, Label, Species, Location, Date,
         Latitude, Longitude, `Accession Number`)  %>% 
  mutate(Latitude = sprintf("%7.4f", Latitude) %>% as.numeric(),
         Longitude = sprintf("%8.4f", Longitude)%>% as.numeric()) %>%
  replace_na(list(Latitude = '-', Longitude = '-')) %>%
  mutate(Nr = row_number()) %>%
  select(Nr, ID:`Accession Number`) %>%
  select(-Label)

export_2_latex(table = table_out,name = 'tables/suppl_tab1.tex')
