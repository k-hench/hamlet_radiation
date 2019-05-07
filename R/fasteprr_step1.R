#!/usr/bin/env Rscript
# Rscript --vanilla fasteprr_step1.R <name.vcf> <OUTPATH> <OUTFILE> <winsize(kb)>
args = commandArgs(trailingOnly=FALSE)
args = args[7:10]
print(args)

library(FastEPRR)

FastEPRR_VCF_step1(vcfFilePath = args[1],
			       erStart = "0.001",
			       qualThreshold = 0,
			       winDXThreshold = 0,
			       winLength = as.character(args[4]),
			       srcOutputFilePath = paste0(args[2],'/',args[3]))
