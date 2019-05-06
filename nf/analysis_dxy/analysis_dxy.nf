#!/usr/bin/env nextflow
// This pipelie includes the anlysis run on the
//   all callable sites data sheet (dxy).

Channel
	.fromFilePairs("../../1_genotyping/3_gatk_filtered/filterd.allBP.vcf.{gz,gz.tbi}")
	.set{ vcf_ch }

Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.set{ lg_ch }

// ------------------------------------
// split the genotypes by LG and reformat the genotypes
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

Channel.from( [[1, "ind"], [2, "may"], [3, "nig"], [4, "pue"], [5, "uni"]] ).into{ bel_spec1_ch; bel_spec2_ch }
Channel.from( [[1, "abe"], [2, "gum"], [3, "nig"], [4, "pue"], [5, "ran"], [6, "uni"]] ).into{ hon_spec1_ch; hon_spec2_ch }
Channel.from( [[1, "nig"], [2, "pue"], [3, "uni"]] ).into{ pan_spec1_ch; pan_spec2_ch }

// Preparation: create all possible species pairs depending on location
//   and combine with genotype subset (for the respective location)

// channel content after joinig: set [0:val(loc), 1:file(vcf), 2:file(pop), 3:val(spec1), 4:val(spec2)]
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

bel_pairs_ch
	.concat( hon_pairs_ch, pan_pairs_ch )
	.combine( geno_ch )
	.set { all_dxy_pairs_ch }

// compute the dxy values along non-overlaping 50kb windows
process dxy_lg {
	label 'L_G32g30h_dxy_lg'
	tag "${spec1}${loc}-${spec2}${loc}_LG${lg}"
	// this process is likely not to finish - somehow the window script
	// fails to finish - I still produces the output though

	input:
	set val( loc ), val( spec1 ), val( spec2 ), val( lg ), file( vcf ), file( geno )  from all_dxy_pairs_ch

	output:
	set val( "${spec1}${loc}-${spec2}${loc}" ), file( "dxy.${spec1}${loc}-${spec2}${loc}.LG${lg}.50kb-5kb.txt.gz" ), val( lg ), val( "${spec1}${loc}" ), val( "${spec2}${loc}" ) into dxy_lg_ch

	script:
	"""
	module load openssl1.0.2
	module load intel17.0.4 intelmpi17.0.4

	vcfsamplenames ${vcf} | \
		awk -v OFS='\t' '{print \$1, substr( \$1, length(\$1) - 5, 6)}' > pop.txt

	mpirun \$NQSII_MPIOPTS -np 1 \
		python \$SFTWR/genomics_general/popgenWindows.py \
		-w 50000 -s 50000 \
		--popsFile pop.txt \
		-p ${spec1}${loc} -p ${spec2}${loc} \
		-g ${geno} \
		-o dxy.${spec1}${loc}-${spec2}${loc}.LG${lg}.50kb-5kb.txt.gz \
		-f phased \
		--writeFailedWindows \
		-T 1
    """
}

dxy_lg_ch
  .groupTuple()
  .set{ tubbled_dxy }

process receive_tuple {
	label 'L_36g47h_receive_tuple'
	publishDir "../../2_analysis/dxy/", mode: 'copy'
	tag "${pop1}-${pop2}"

	input:
	set comp, dxy, lg, pop1, pop2 from tubbled_dxy

	output:
	file( "dxy.${pop1}-${pop2}.50kb-5kb.tsv.gz" ) into dxy_output_ch

	script:
	"""
	zcat dxy.${pop1}-${pop2}.LG01.50kb-5kb.txt.gz | \
	head -n 1 > dxy.${pop1}-${pop2}.50kb-5kb.tsv;

	for j in {01..24};do
		echo "-> LG\$j"
		zcat dxy.${pop1}-${pop2}.LG\$j.50kb-5kb.txt.gz | \
			awk 'NR>1{print}' >> dxy.${pop1}-${pop2}.50kb-5kb.tsv;
	done

	gzip dxy.${pop1}-${pop2}.50kb-5kb.tsv
	"""
}
