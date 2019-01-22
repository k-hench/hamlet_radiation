#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla assign_phenotypes.R ${fam_file} ${trait}
# ===============================================================
# This script 
# ---------------------------------------------------------------
# ===============================================================
# args <- c('plink_test_binary.fam', 'inds.txt', 'pheno1')
args = commandArgs(trailingOnly=FALSE)
args = args[7:9]
print(args)

infile <- as.character(args[1])
pheno_file <- as.character(args[2])
trait <- as.character(args[3])
library(tidyverse)

phenotypes <- read_delim(file = pheno_file,delim = '\t')


data <- read_delim(infile,delim = ' ',
                   col_names = c("ID", "Within_family_ID", "ID_father", "ID_mother", "Sex", "Phenotype")) %>%
  left_join(phenotypes) %>%
  select(ID:Sex,trait)

write_delim(x = data, path = infile, delim = ' ', col_names = FALSE)
