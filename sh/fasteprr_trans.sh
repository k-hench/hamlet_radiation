#/usr/bin/bash
# fastprr_trans.sh <infile> <lg> <outfile>

cat $1 | \
    paste - - | \
    sed 's/ /\t/g; s/:\tRho:/\t/g; s/-/\t/g' | \
    cut -f 2,3,4 | \
    awk -v OFS="\t" -v lg=$2 'BEGIN{print "CHROM\tBIN_START\tBIN_END\tRHO"}NR>0{print lg,$1*1000,$2*1000,$3}' > $3.rho.txt
