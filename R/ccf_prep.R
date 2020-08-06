#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/ccf_prep.R 2_analysis/geva/ LG01 16_21-30nigpan 20553puehon
# ===============================================================
# This script 
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "2_analysis/geva/", "LG01", "16_21-30nigpan", "20553puehon")
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(tidyverse)
library(vroom)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

cli::rule( left = str_c(crayon::bold('Script: '),crayon::red(script_name)))
args = args[7:length(args)]
cat(' ')
cat(str_c(crayon::green(cli::symbol$star),' ', 1:length(args),': ',crayon::green(args),'\n'))
cli::rule(right = getwd())

# config -----------------------
geva_path <- as.character(args[1])
chrom <- as.character(args[2])
target <- as.character(args[3])
querry <- as.character(args[4])

data <- vroom::vroom("~/work/tests/ancient_allele_vcf/phased2.vcf", delim = "\t",skip = 5)

data_geva <-   vroom::vroom(file = str_c(geva_path, chrom,".sites.txt.gz"), delim = " ") %>%
  left_join(vroom::vroom(str_c(geva_path, chrom,".marker.txt.gz"), delim = " ")) %>%
  filter(Clock == "J", Filtered == 0) %>%
  arrange(Position) %>%
  mutate(`#CHROM` = str_c("LG", str_pad(Chromosome,width = 2,pad = "0")),
         POS = Position) %>%
  select(`#CHROM`,POS,PostMedian)

split_id <- function(data, target, querry){
  target <- enquo(target)
  querry <- enquo(querry)
  
  data %>%
    separate(!!target, into = c("target_A", "target_B"), sep = "\\|",convert = TRUE) %>%
    separate(!!querry, into = c("querry_A", "querry_B"), sep = "\\|",convert = TRUE) %>%
    mutate(tAtB = as.numeric((target_A-target_B) == 0),
           tAqA = as.numeric((target_A-querry_A) == 0),
           tAqB = as.numeric((target_A-querry_B) == 0),
           tBqA = as.numeric((target_B-querry_A) == 0),
           tBqB = as.numeric((target_B-querry_B) == 0),
           qAqB = as.numeric((querry_A-querry_B) == 0)) %>%
    select(`#CHROM`, POS, tAtB:qAqB) %>%
    left_join(data_geva, by = c(`#CHROM` = "#CHROM", POS = "POS")) %>%
    filter(!is.na(PostMedian)) %>%
    arrange(PostMedian) 
}

data %>% 
  split_id(target, querry) %>%
    set_names(nm = c("#CHROM", "POS",
                     str_c(target, "_1-", target, "_2"),
                     str_c(target, "_1-", querry, "_1"),
                     str_c(target, "_1-", querry, "_2"),
                     str_c(target, "_2-", querry, "_1"),
                     str_c(target, "_2-", querry, "_2"),
                     str_c(querry, "_1-", querry, "_2"),
                     "PostMedian")) %>%
  write_tsv(str_c(chrom,"_",target, "-", querry,"_ccf_prep.tsv"))

system(str_c("gzip ",chrom,"_",target, "-", querry,"_ccf_prep.tsv"))

