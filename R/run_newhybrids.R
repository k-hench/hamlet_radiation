#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla run_newhybrids.R
# ===============================================================
# This script
#   ---------------------------------------------------------------
# ===============================================================
# args <- c()
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(parallelnewhybrid)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

args <- process_input(script_name, args)
# config -----------------------

## Get the file path to the working directory, will be used to allow a universal example
parallelnh_LINUX(folder.data = str_c(getwd(),"/nh_input/"),
                 where.NH = "/gss/work/foge1444/software/newhybrids/",
                 burnin = 1000000,
                 sweeps = 10000000)
