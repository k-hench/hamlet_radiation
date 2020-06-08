#!/usr/bin/env nextflow
// git 9.1
// open genotype data
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.into{ vcf_loc1; vcf_loc2; vcf_loc3 }

// git 9.2
// initialize location channel
Channel
	.from( "bel", "hon", "pan")
	.set{ locations_ch }

// git 9.3
// define location specific sepcies set
Channel.from( [[1, "ind"], [2, "may"], [3, "nig"], [4, "pue"], [5, "uni"]] ).into{ bel_spec1_ch; bel_spec2_ch }
Channel.from( [[1, "abe"], [2, "gum"], [3, "nig"], [4, "pue"], [5, "ran"], [6, "uni"]] ).into{ hon_spec1_ch; hon_spec2_ch }
Channel.from( [[1, "nig"], [2, "pue"], [3, "uni"]] ).into{ pan_spec1_ch; pan_spec2_ch }

// git 9.4
// prepare pairwise new_hybrids
// ------------------------------
/* (create all possible species pairs depending on location
   and combine with genotype subset (for the respective location))*/
bel_pairs_ch = Channel.from( "bel" )
    .combine( vcf_loc1 )
    .combine(bel_spec1_ch)
    .combine(bel_spec2_ch)
    .filter{ it[3] < it[5] }
    .map{ it[0,1,2,4,6]}
hon_pairs_ch = Channel.from( "hon" )
    .combine( vcf_loc2 )
    .combine(hon_spec1_ch)
    .combine(hon_spec2_ch)
    .filter{ it[3] < it[5] }
    .map{ it[0,1,2,4,6]}
pan_pairs_ch = Channel.from( "pan" )
    .combine( vcf_loc3 )
    .combine(pan_spec1_ch)
    .combine(pan_spec2_ch)
    .filter{ it[3] < it[5] }
    .map{ it[0,1,2,4,6]}
bel_pairs_ch.concat( hon_pairs_ch, pan_pairs_ch  ).set { all_fst_pairs_ch }

// git 9.5
// comute pairwise fsts
process fst_run {
	label 'L_20g45m_fst_run'
	tag "${spec1}${loc}-${spec2}${loc}"

	input:
	set val( loc ), val( vcfidx ), file( vcf ), val( spec1 ), val( spec2 ) from all_fst_pairs_ch

	output:
	set val( loc ), val( spec1 ), val( spec2 ), file( "${vcf[0]}" ), file( "*.fst.tsv.gz" ), file( "${spec1}${loc}.pop"), file( "${spec2}${loc}.pop") into fst_SNPS

	script:
	"""
	vcfsamplenames ${vcf[0]} | grep ${spec1}${loc} > ${spec1}${loc}.pop
	vcfsamplenames ${vcf[0]} | grep ${spec2}${loc} > ${spec2}${loc}.pop

	vcftools --gzvcf ${vcf[0]} \
		 --weir-fst-pop ${spec1}${loc}.pop \
		 --weir-fst-pop ${spec2}${loc}.pop \
		 --stdout | gzip > ${spec1}${loc}-${spec2}${loc}.fst.tsv.gz
	"""
}

// git 9.6
process filter_fst {
	label 'L_8g15m_filter_fst'
	tag "${spec1}${loc}-${spec2}${loc}"

	input:
	set val( loc ), val( spec1 ), val( spec2 ), file( vcf ), file( fst ), file( pop1 ), file( pop2 ) from fst_SNPS

	output:
	set val( loc ), val( spec1 ), val( spec2 ), file( vcf ), file( pop1 ), file( pop2 ), file( "*SNPs.snps" ) into filter_SNPs


	script:
	"""
	Rscript --vanilla \$BASE_DIR/R/filter_snps.R ${fst} 800 ${spec1}${loc}-${spec2}${loc}
	"""
}

// git 9.7
process prep_nh_input {
	label 'L_8g15m_prep_nh'
	tag "${spec1}${loc}-${spec2}${loc}"

	input:
	set val( loc ), val( spec1 ), val( spec2 ), file( vcf ), file( pop1 ), file( pop2 ), file( snps ) from filter_SNPs


	output:
	set val( loc ), val( spec1 ), val( spec2 ), file( "*_individuals.txt" ), file( "*.80SNPs.txt")  into newhybrids_input


	script:
	"""
	vcftools \
  --gzvcf ${vcf} \
	--keep ${pop1} \
	--keep ${pop2} \
	--thin 5000 \
	--out newHyb.${spec1}${loc}-${spec2}${loc} \
	--positions ${snps} \
	--recode

	grep '#' newHyb.${spec1}${loc}-${spec2}${loc}.recode.vcf > newHyb.${spec1}${loc}-${spec2}${loc}.80SNPs.vcf
	grep -v '#' newHyb.${spec1}${loc}-${spec2}${loc}.recode.vcf | \
		shuf -n 80 | \
		sort -k 1 -k2 >> newHyb.${spec1}${loc}-${spec2}${loc}.80SNPs.vcf

	grep '#CHROM' newHyb.${spec1}${loc}-${spec2}${loc}.80SNPs.vcf | \
		cut -f 10- | \
		sed 's/\\t/\\n/g' > newHyb.${spec1}${loc}-${spec2}${loc}.80SNPs_individuals.txt

	/usr/bin/java -Xmx1024m -Xms512M \
		-jar \$SFTWR/PGDSpider/PGDSpider2-cli.jar \
		-inputfile newHyb.${spec1}${loc}-${spec2}${loc}.80SNPs.vcf \
		-inputformat VCF \
		-outputfile newHyb.${spec1}${loc}-${spec2}${loc}.80SNPs.txt \
		-outputformat NEWHYBRIDS \
		-spid \$BASE_DIR/ressources/vcf2nh.spid
	"""
}

// git 9.8
// copy of nh_input is needed because nh can't read links...
process run_nh {
	label 'L_20g15h4x_run_nh'
	tag "${spec1}${loc}-${spec2}${loc}"
	publishDir "../../2_analysis/newhyb/", mode: 'copy'

	input:
	set val( loc ), val( spec1 ), val( spec2 ), file( inds ), file( snps ) from newhybrids_input

	output:
	set file( "nh_input/NH.Results/newHyb.*/*_individuals.txt" ), file( "nh_input/NH.Results/newHyb.*/*_PofZ.txt" )  into newhybrids_output

	script:
	"""
	mkdir -p nh_input
	cp ${snps} nh_input/${snps}
	cp ${inds} nh_input/${inds}

	Rscript --vanilla \$BASE_DIR/R/run_newhybrids.R
	"""
}
