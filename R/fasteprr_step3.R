#!/usr/bin/env Rscript
# Rscript --vanilla fasteprr_step3.R <INPATH> <OUTPATH2> <OUTPATH3>
args = commandArgs(trailingOnly=FALSE)
args = args[7:9]

print(args)
library(FastEPRR)

FastEPRR_VCF_step3(srcFolderPath=args[1], DXFolderPath=args[2], finalOutputFolderPath=args[3])  

