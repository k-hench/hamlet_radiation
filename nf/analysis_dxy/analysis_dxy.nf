#!/usr/bin/env nextflow
// This pipelie includes the anlysis run on the
//   all callable sites data sheet (dxy).

// git 4.1
// load genotypes
Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/filterd.allBP.vcf.{gz,gz.tbi}")
	.into{ vcf_ch; vcf_pi_ch }

// git 4.2
// initialize LGs
Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.set{ lg_ch }

// git 4.3
// split by LG and reformat the genotypes
process split_allBP {
	label 'L_32g15h_split_allBP'
	tag "LG${lg}"

	input:
	set val( lg ), vcfId, file( vcf ) from lg_ch.combine( vcf_ch )

	output:
	set val( lg ), file( 'filterd.allBP.vcf.gz' ), file( "allBP.LG${lg}.geno.gz" ) into geno_ch

	script:
	"""
	module load openssl1.0.2

	vcftools --gzvcf ${vcf[0]} \
		--chr LG${lg} \
		--recode \
		--stdout | bgzip  > allBP.LG${lg}.vcf.gz

	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py \
		-i allBP.LG${lg}.vcf.gz  | gzip > allBP.LG${lg}.geno.gz
	"""
}

// git 4.4
// define location specific sepcies set
Channel.from( [[1, "ind"], [2, "may"], [3, "nig"], [4, "pue"], [5, "uni"]] ).into{ bel_spec1_ch; bel_spec2_ch }
Channel.from( [[1, "abe"], [2, "gum"], [3, "nig"], [4, "pue"], [5, "ran"], [6, "uni"]] ).into{ hon_spec1_ch; hon_spec2_ch }
Channel.from( [[1, "nig"], [2, "pue"], [3, "uni"]] ).into{ pan_spec1_ch; pan_spec2_ch }

// git 4.5
// init all sampled populations (for pi)
Channel
	.from('indbel', 'maybel', 'nigbel', 'puebel', 'unibel', 'abehon', 'gumhon', 'nighon', 'puehon', 'ranhon', 'unihon', 'nigpan', 'puepan', 'unipan')
	.set{spec_dxy}


// git 4.6
Channel
	.from( 1, 5 )
	.into{ kb_ch; kb_ch2; kb_ch3 }

// git 4.7
// prepare pair wise dxy
// ------------------------------
// create all possible species pairs depending on location
//   and combine with genotype subset (for the respective location)
// ------------------------------
// channel content after joinig:
// set [0:val(loc), 1:file(vcf), 2:file(pop), 3:val(spec1), 4:val(spec2)]
// ------------------------------
bel_pairs_ch = Channel.from( "bel" )
	.combine( bel_spec1_ch )
	.combine( bel_spec2_ch )
	.filter{ it[1] < it[3] }
	.map{ it[0,2,4]}
hon_pairs_ch = Channel.from( "hon" )
	.combine( hon_spec1_ch )
	.combine(hon_spec2_ch)
	.filter{ it[1] < it[3] }
	.map{ it[0,2,4]}
pan_pairs_ch = Channel.from( "pan" )
	.combine( pan_spec1_ch )
	.combine(pan_spec2_ch)
	.filter{ it[1] < it[3] }
	.map{ it[0,2,4]}

// git 4.8
// combine species pair with genotypes (and window size)
bel_pairs_ch
	.concat( hon_pairs_ch, pan_pairs_ch )
	.combine( geno_ch )
	.combine( kb_ch )
	.into { all_dxy_pairs_ch; random_dxy_pairs_ch }

// git 4.9
// compute the dxy values
process dxy_lg {
	label 'L_G32g15h_dxy_lg'
	tag "${spec1}${loc}-${spec2}${loc}_LG${lg}"

	input:
	set val( loc ), val( spec1 ), val( spec2 ), val( lg ), file( vcf ), file( geno ), val( kb ) from all_dxy_pairs_ch

	output:
	set val( "${spec1}${loc}-${spec2}${loc}-${kb}" ), file( "dxy.${spec1}${loc}-${spec2}${loc}.LG${lg}.${kb}0kb-${kb}kb.txt.gz" ), val( lg ), val( "${spec1}${loc}" ), val( "${spec2}${loc}" ), val( kb ) into dxy_lg_ch

	script:
	"""
	module load openssl1.0.2
	module load intel17.0.4 intelmpi17.0.4

	zcat ${geno} | \
		head -n 1 | \
		cut -f 3- | \
		sed 's/\\t/\\n/g' | \
		awk -v OFS='\\t' '{print \$1, substr( \$1, length(\$1) - 5, 6)}' > pop.txt

	mpirun \$NQSII_MPIOPTS -np 1 \
		python \$SFTWR/genomics_general/popgenWindows.py \
		-w ${kb}0000 -s ${kb}000 \
		--popsFile pop.txt \
		-p ${spec1}${loc} -p ${spec2}${loc} \
		-g ${geno} \
		-o dxy.${spec1}${loc}-${spec2}${loc}.LG${lg}.${kb}0kb-${kb}kb.txt.gz \
		-f phased \
		--writeFailedWindows \
		-T 1
    """
}

// git 4.10
// collect all LGs for each species pair
dxy_lg_ch
  .groupTuple()
  .set{ tubbled_dxy }

// git 4.11
// concatinate all LGs for each species pair
process receive_tuple {
	label 'L_20g2h_receive_tuple'
	publishDir "../../2_analysis/dxy/${kb[0]}0k/", mode: 'copy'
	tag "${pop1[0]}-${pop2[0]}"

	input:
	set val( comp ), file( dxy ), val( lg ), val( pop1 ), val( pop2 ), val( kb ) from tubbled_dxy

	output:
	file( "dxy.${pop1[0]}-${pop2[0]}.${kb[0]}0kb-${kb[0]}kb.tsv.gz" ) into dxy_output_ch

	script:
	"""
	zcat dxy.${pop1[0]}-${pop2[0]}.LG01.${kb[0]}0kb-${kb[0]}kb.txt.gz | \
	head -n 1 > dxy.${pop1[0]}-${pop2[0]}.${kb[0]}0kb-${kb[0]}kb.tsv;

	for j in {01..24};do
		echo "-> LG\$j"
		zcat dxy.${pop1[0]}-${pop2[0]}.LG\$j.${kb[0]}0kb-${kb[0]}kb.txt.gz | \
			awk 'NR>1{print}' >> dxy.${pop1[0]}-${pop2[0]}.${kb[0]}0kb-${kb[0]}kb.tsv;
	done

	gzip dxy.${pop1[0]}-${pop2[0]}.${kb[0]}0kb-${kb[0]}kb.tsv
	"""
}

// git 4.12
// collect a species pair to randomize
Channel
	.from( [['bel', 'ind', 'may']] )
	.set{ random_run_ch }

// git 4.13
// setup channel contnent for random channel
Channel
	.from( 1 )
	.combine( random_run_ch )
	.combine( kb_ch2 )
	.filter{ it[4] == 5 }
	.set{ random_sets_ch }

// git 4.14
// permute the population assignment (the randomization)
process randomize_samples {
	label 'L_20g15h_randomize_samples'
	publishDir "../../2_analysis/fst/${kb}0k/random", mode: 'copy' , pattern: "*_windowed.weir.fst.gz"
	module "R3.5.2"

	input:
	set val( random_set ), val( loc ), val(spec1), val(spec2), val( kb ) from random_sets_ch

	output:
	set random_set, file( "random_pop.txt" ) into random_pops_ch
	file( "*_windowed.weir.fst.gz") into random_fst_out

	script:
	"""
	cut -f 2,3 \$BASE_DIR/metadata/sample_info.txt | \
		grep "${loc}" | \
		grep "${spec1}\\|${spec2}" > pop_prep.tsv

	Rscript --vanilla \$BASE_DIR/R/randomize_pops.R

	grep A random_pop.txt | cut -f 1  > pop1.txt
	grep B random_pop.txt | cut -f 1  > pop2.txt

	vcftools \
	  --gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered/filterd_bi-allelic.allBP.vcf.gz \
	  --weir-fst-pop pop1.txt \
	  --weir-fst-pop pop2.txt \
	  --fst-window-step ${kb}0000 \
	  --fst-window-size ${kb}0000 \
	  --stdout | gzip > ${loc}-aaa-bbb.${kb}0k.random_${spec1}_${spec2}_windowed.weir.fst.gz

	"""
}

// git 4.15
// pick random pair of interest
random_dxy_pairs_ch
	.filter{ it[0] == 'bel' && it[1] == 'ind' && it[2] == 'may'  && it[6] == 5 }
	.combine( random_pops_ch )
	.set{ random_assigned_ch }

// git 4.16
// compute the dxy values
process dxy_lg_random {
	label 'L_G32g15h_dxy_lg_random'
	tag "aaa${loc}-bbb${loc}_LG${lg}"
	module "R3.5.2"

	input:
	set val( loc ), val( spec1 ), val( spec2 ), val( lg ), file( vcf ), file( geno ), val( kb ), val( random_set ), file( pop_file ) from random_assigned_ch

	output:
	set val( "aaa${loc}-bbb${loc}-${kb}0kb" ), file( "dxy.aaa${loc}-bbb${loc}.LG${lg}.${kb}0kb-${kb}kb.txt.gz" ), val( lg ), val( "aaa${loc}" ), val( "bbb${loc}" ), val( kb ) into dxy_random_lg_ch

	script:
	"""
	module load openssl1.0.2
	module load intel17.0.4 intelmpi17.0.4

	mpirun \$NQSII_MPIOPTS -np 1 \
		python \$SFTWR/genomics_general/popgenWindows.py \
		-w ${kb}0000 -s ${kb}000 \
		--popsFile ${pop_file} \
		-p A -p B \
		-g ${geno} \
		-o dxy.aaa${loc}-bbb${loc}.LG${lg}.${kb}0kb-${kb}kb.txt.gz \
		-f phased \
		--writeFailedWindows \
		-T 1
	 """
}

// git 4.17
// collect all LGs of random run
dxy_random_lg_ch
.groupTuple()
.set{ tubbled_random_dxy }

// git 4.18
// concatinate all LGs of random run
process receive_random_tuple {
	label 'L_20g2h_receive_random_tuple'
	publishDir "../../2_analysis/dxy/random/", mode: 'copy'

	input:
	set val( comp ), file( dxy ), val( lg ), val( pop1 ), val( pop2 ), val( kb )  from tubbled_random_dxy

	output:
	file( "dxy.${pop1[0]}-${pop2[0]}.${kb[0]}0kb-${kb[0]}kb.tsv.gz" ) into dxy_random_output_ch

	script:
	"""
	zcat dxy.${pop1[0]}-${pop2[0]}.LG01.${kb[0]}0kb-${kb[0]}kb.txt.gz | \
	head -n 1 > dxy.${pop1[0]}-${pop2[0]}.${kb[0]}0kb-${kb[0]}kb.tsv;

	for j in {01..24};do
		echo "-> LG\$j"
		zcat dxy.${pop1[0]}-${pop2[0]}.LG\$j.${kb[0]}0kb-${kb[0]}kb.txt.gz | \
			awk 'NR>1{print}' >> dxy.${pop1[0]}-${pop2[0]}.${kb[0]}0kb-${kb[0]}kb.tsv;
	done

	gzip dxy.${pop1[0]}-${pop2[0]}.${kb[0]}0kb-${kb[0]}kb.tsv
	"""
}

// ---------------------------------------------------------------
// The pi part need to be run AFTER the global fst outlier
// windows were selected (REMEMBER TO CHECK FST OUTLIER DIRECTORY)
// ---------------------------------------------------------------

// git 4.19
// calculate pi per species
process pi_per_spec {
	label 'L_32g15h_pi'
	tag "${spec}"
	publishDir "../../2_analysis/pi/${kb}0k", mode: 'copy'

	input:
	set val( spec ), vcfId, file( vcf ), val( kb ) from spec_dxy.combine( vcf_pi_ch ).combine( kb_ch3 )

	output:
   file( "*.${kb}0k.windowed.pi.gz" ) into pi_50k

	script:
	"""
	module load openssl1.0.2

	vcfsamplenames ${vcf[0]} | \
		grep ${spec} > pop.txt

	vcftools --gzvcf ${vcf[0]} \
		--keep pop.txt \
		--window-pi ${kb}0000 \
		--window-pi-step ${kb}000 \
		--out ${spec}.${kb}0k 2> ${spec}.pi.log
	gzip ${spec}.${kb}0k.windowed.pi

	tail -n +2 \$BASE_DIR/2_analysis/summaries/fst_outliers_998.tsv | \
		cut -f 2,3,4 > outlier.bed

	vcftools --gzvcf ${vcf[0]} \
		--keep pop.txt \
		--exclude-bed outlier.bed \
		--window-pi ${kb}0000 \
		--window-pi-step ${kb}000\
		--out ${spec}_no_outlier.${kb}0k 2> ${spec}_${kb}0k_no_outllier.pi.log
	gzip ${spec}_no_outlier.${kb}0k.windowed.pi
	"""
}