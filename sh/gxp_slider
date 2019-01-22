#!/usr/bin/bash
# by: Kosmas Hench: 2019-01-22
# usage gxp_slider <in.GxP.txt.gz> <Windowsize (bp)> <Increment (bp)>
# -------------------------------------------------------------------------------------------------------------------
# !! beware of awk array size maximum - if the sequence lenght is to big or if the window size & step are to small !!
# !!   (hence, if too many windows are created, the array index maximum is early entries are being overwritten     !!
# -------------------------------------------------------------------------------------------------------------------
# $1 = infile; $2 = windowsize; $3 = increment
j=$(echo $1 | sed 's/.GxP.txt.gz//g')
WIN=$2;
STEP=$3;
WIN_LAB=$(($WIN/1000));
STEP_LAB=$(($STEP/1000));
RED='\033[0;31m';
NC='\033[0m';

# print the header (the first line of input)
# and then run the specified command on the body (the rest of the input)
# use it in a pipeline, e.g. ps | body grep somepattern
body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}

for k in {01..24};do
  echo -e  "--- ${RED}$j${NC} - $k ---"

	zcat $1 | head -n 1  | cut -f 1,2,6-8 > $j.LG$k.tmp

	zcat $1 | \
		awk -v k=LG$k -v OFS='\t' '$1==k {print}' | \
		cut -f 1,2,6-8 | \
		awk -v OFS='\t' '{$3=-log($3)/log(10);$4=-log($4)/log(10);$5=-log($5)/log(10);print $0}' >> $j.LG$k.tmp

	grep -v nan $j.LG$k.tmp | \
	awk -v OFS="\t" -v w=$WIN -v s=$STEP -v r=$k 'BEGIN{window=w;slide=s;g=0;OL=int(window/slide);}
	{if(NR==1 && r==01){print "CHROM","BIN_START","BIN_END","N_SNPs","MID_POS","BIN_RANK","BIN_NR","SNP_DENSITY","AVG_"$3 ,"AVG_"$4,"AVG_"$5;}}
	{if(NR>1){g=int(($2-1)/slide)+1;{for (i=0;i<=OL-1;i++){if(g-i >0){A[g-i]+=$2; B[g-i]++;C[g-i]+=$3;D[g-i]+=$4;E[g-i]+=$5;G[g-i]=g-i;H[g-i]=$1;}}}};}
	END{for(i in A){print H[i],(G[i]-1)*slide,(G[i]-1)*slide+window,B[i],A[i]/B[i],G[i],k+1,B[i]/window,C[i]/B[i],D[i]/B[i],E[i]/B[i];k++}}' |
	body sort -k2 -n | \
	sed -r '/^\s*$/d' >> ${j}.${WIN_LAB}k.${STEP_LAB}k.txt
done

gzip ${j}.${WIN_LAB}k.${STEP_LAB}k.txt
rm ./$j*.tmp
