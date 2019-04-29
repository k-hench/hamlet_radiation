#!/usr/bin/env bash
BASE_NAME=$(echo $1 | sed 's/.vcf.gz//g')

zcat $1 | \
grep -v '^#' | \
cut -f 1,2  | \
awk -v nsps=$2 -v OFS="\t" \
'BEGIN{print "lg","start","end","n_snps","win_id","gwin_id"}
{if(NR == 1){lg=$1; start = $2; end = $2; snpcount = 1; wincount = 1; gwincount = 1}
else if(lg != $1){print lg, start, end, snpcount, wincount, gwincount; lg = $1; start = $2; end = $2; snpcount = 1; wincount = 1; gwincount = ++gwincount}
else if(lg == $1){
	  if(snpcount < nsps){end = $2; snpcount = ++snpcount}
  else if(snpcount == nsps){print lg, start, end, snpcount, wincount, gwincount; start = $2; end = $2; snpcount = 1; wincount = ++wincount; gwincount = ++gwincount}
  else {print "ERROR: SOMETHING WENT WRONG!!" }
  }}
  END {print lg, start, end, snpcount, wincount, gwincount}' > .tmp.tsv

awk '{if($2 != $3){print}}' .tmp.tsv > $BASE_NAME.snp_windows.$2.tsv
rm .tmp.tsv