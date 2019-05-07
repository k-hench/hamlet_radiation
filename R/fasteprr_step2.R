#!/usr/bin/env Rscript
# Rscript --vanilla fasteprr_step2.R <INPATH> <OUTPATH> <JOBSTART>
args = commandArgs(trailingOnly=FALSE)
args = args[7:9]

print(args)
library(FastEPRR)

j <- as.numeric(args[3])
FastEPRR_VCF_step2(srcFolderPath=args[1], jobNumber=250, currJob=j, DXOutputFolderPath=args[2])

