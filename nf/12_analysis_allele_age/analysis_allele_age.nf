#!/usr/bin/env nextflow
// git 12.1
// open genotype data
Channel
	.fromFilePairs("../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}")
	.set{ vcf_ch }

// git 12.2
// initialize LGs
Channel
	.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
	.map{"LG" + it}
	.combine( vcf_ch )
	.set{ lg_ch }

// git 12.3
// subset the genotypes by location
process prepare_vcf {
	label "L_20g2h_prepare_vcf"

	input:
	set val( lg ), val( vcfidx ), file( vcf ) from lg_ch

	output:
	set val( lg ), file( "${lg}_integer.vcf" )  into ( vcf_prep_ch )

	script:
	"""
	# subset by LG and add AC info field
	vcftools \
	  --gzvcf ${vcf[0]} \
		--chr ${lg} \
		--recode \
		--stdout | \
		sed 's/\\(##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\\)/\\1\\n##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">/' | \
		bgzip > ${lg}.vcf.gz

	# determine ancestral state based on invariant sites in outgoup
	zcat ${lg}.vcf.gz | \
		grep -v "^#" | \
		awk -v OFS="\\t" \
		'BEGIN{print "#CHROM","FROM","TO","AA"}
		{o1=substr(\$107, 1, 1);
		o2=substr(\$107, 3, 3);
		o3=substr(\$167, 1, 1);
		o4=substr(\$167, 3, 3);
		o5=substr(\$179, 1, 1);
		o6=substr(\$179, 3, 3);
		if (o1 == o2 && o3 == o4 && o5 == o6 && o1 == o3 && o1 == o5){
		aa = \$(4+o1)} else {aa = "."};
		print \$1,\$2,\$2,aa}' > ${lg}_annotations.bed

	# determine allele frquencies
	vcftools \
		--gzvcf ${lg}.vcf.gz \
		--freq \
		--stdout | \
		sed 's/{ALLELE:FREQ}/ALLELE1\\tALLELE2/' > ${lg}_allele_counts.tsv

	# determine ancestral state for variant sites in outgoup based on allele freq
	Rscript --vanilla \$BASE_DIR/R/major_allele.R ${lg}_allele_counts.tsv ${lg}_annotations.bed

	bgzip ${lg}_annotations_maj.bed

	tabix -s 1 -b 2 -e 3 ${lg}_annotations_maj.bed.gz

	# add ancestral state annotation
	zcat ${lg}.vcf.gz | \
		vcf-annotate -a ${lg}_annotations_maj.bed.gz \
		-d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
		-c CHROM,FROM,TO,INFO/AA | \
		sed 's/LG//g'  \
		> ${lg}_integer.vcf
	"""
}

// git 12.4
// set ancestral state in vcf
process set_ancestral_states {
	label 'L_2g15m_ancestral_states'
	publishDir "../../1_genotyping/5_ancestral_allele", mode: 'copy'

	input:
	set val( lg ), file( vcf ) from ( vcf_prep_ch )

	output:
	set val( lg ), file( "${lg}_aa.vcf.gz" ) into ( vcf_aa_ch )

	script:
	"""
	java -jar \$SFTWR/jvarkit/dist/vcffilterjdk.jar \
		-f \$BASE_DIR/js/script.js ${vcf} | \
		bgzip > ${lg}_aa.vcf.gz
	"""
}

// git 12.5
// set ancestral state in vcf
process create_positions {
	label 'L_20g2h_create_positions'
	publishDir "../../2_analysis/sliding_phylo/", mode: 'copy'

	input:
	set val( lg ), file( vcf ) from ( vcf_aa_ch )

	output:
	set val( lg ), file( "${lg}_aa_h_variant.vcf.gz" ), file( "${lg}_positions.txt" ) into ( positions_ch )

	script:
	"""
	echo -e "20478tabhon\\n28393torpan\\ns_tort_3torpan" > outgr.pop

	# keeping only sites that are variant within hamlets
	vcftools \
		--gzvcf ${vcf} \
		--remove outgr.pop \
		--recode \
		--stdout | \
		vcftools \
		--gzvcf - \
		--mac 1 \
		--recode \
		--stdout | \
		bgzip > ${lg}_aa_no_outgroup.vcf.gz

	zcat ${lg}_aa_no_outgroup.vcf.gz | \
		grep -v "^#" | \
		cut -f 1,2 | \
		head -n -1 > ${lg}_positions_prep.txt

	vcftools \
		--gzvcf ${vcf} \
		--positions ${lg}_positions_prep.txt \
		--recode \
		--stdout | \
		bgzip > ${lg}_aa_h_variant.vcf.gz

	cut -f 2 ${lg}_positions_prep.txt > ${lg}_positions.txt
	"""
}


// git 12.5
// set ancestral state in vcf
// all in one takes too long...
process pre_split {
	label 'L_2g2h_pre_split'
	publishDir "../../2_analysis/geva/", mode: 'copy'

	input:
	set val( lg ), file( vcf ), file( pos ) from ( positions_ch )

	output:
	set val( lg ), file( "pre_positions/pre_*" ), file( "*.bin" ), file( "*.marker.txt" ), file( "*.sample.txt" ) into ( geva_setup_ch )
	set val( lg ), file( "inner_pos.txt" ), file( vcf ) into ( ccf_vcf_ch )

	script:
	"""
	mkdir -p pre_positions

	head -n -1 ${pos} | \
	 tail -n +2  > inner_pos.txt

	split inner_pos.txt -a 4 -l 25000 -d pre_positions/pre_

	r=\$(awk -v k=${lg} '\$1 == k {print \$4}' \$BASE_DIR/ressources/avg_rho_by_LG.tsv)

	geva_v1beta \
		--vcf ${vcf} --rec \$r --out ${lg}
	"""
}

// git 12.6
// set ancestral state in vcf
process run_geva {
	label 'L_30g15h6x_run_geva'

	input:
	set val( lg ), file( pos ), file( bin ), file( marker ), file( sample ) from geva_setup_ch.transpose()

	output:
	set val( lg ), file( "*.sites.txt.gz" ), file( "*.pairs.txt.gz" ) into ( output_split_ch )

	script:
	"""
  pref=\$(echo "${pos}" | sed 's=^.*/A==; s=pre_positions/pre_==')

	mkdir -p sub_positions sub_results

	split ${pos} -a 4 -l 250 -d sub_positions/sub_pos_\${pref}_

	r=\$(awk -v k=${lg} '\$1 == k {print \$4}' \$BASE_DIR/ressources/avg_rho_by_LG.tsv)

	for sp in \$(ls sub_positions/sub_pos_\${pref}_*); do
		run_id=\$(echo \$sp | sed "s=sub_positions/sub_pos_\${pref}_==")

		geva_v1beta \
			 -t 6 \
			 -i ${bin} \
			 -o sub_results/${lg}_\${pref}_\${run_id}\
			 --positions \$sp \
			 --Ne 30000 \
			 --mut 3.7e-08 \
			 --hmm \$SFTWR/geva/hmm/hmm_initial_probs.txt \$SFTWR/geva/hmm/hmm_emission_probs.txt

		tail -n +2 sub_results/${lg}_\${pref}_\${run_id}.sites.txt >> ${lg}_\${pref}.sites.txt
		tail -n +2 sub_results/${lg}_\${pref}_\${run_id}.pairs.txt >> ${lg}_\${pref}.pairs.txt
	done

	gzip ${lg}_\${pref}.sites.txt
	gzip ${lg}_\${pref}.pairs.txt
	"""
}

// git 12.7
// collect by lg
process collect_by_lg {
	label 'L_2g2h_collect'
	publishDir "../../2_analysis/geva/", mode: 'copy'

	input:
	set val( lg ), file( sites ), file( pairs ) from output_split_ch.groupTuple()

	output:
	set val( lg ), file( "*.sites.txt.gz" ), file( "*.pairs.txt.gz" ) into ( output_lg_ch )

	script:
	"""
	echo "MarkerID Clock Filtered N_Concordant N_Discordant PostMean PostMode PostMedian" > ${lg}.sites.txt
	echo "MarkerID Clock SampleID0 Chr0 SampleID1 Chr1 Shared Pass SegmentLHS SegmentRHS Shape Rate" > ${lg}.pairs.txt

	zcat ${sites}  >> ${lg}.sites.txt
	zcat ${pairs}  >> ${lg}.pairs.txt

	gzip ${lg}.sites.txt
	gzip ${lg}.pairs.txt
	"""
}

// git 12.8
// load ccf pairs
Channel
	.fromPath('../../ressources/plugin/ccf_pairs.tsv')
	.splitCsv(header:true, sep:"\t")
	.map{ row -> [ target:row.target, querry:row.querry] }
	.combine( output_lg_ch.join(ccf_vcf_ch) )
	.set { ccf_pairs_ch }

// git 12.9
// prepare ccf pair
process prep_ccf_pair {
	label 'L_2g2h_prep_ccf_pair'
	publishDir "../../2_analysis/ccf/", mode: 'copy'

	input:
	set val( ccfpair ), val( lg ), file( sites ), file( pairs ), file( pos ), file( vcf ) from ccf_pairs_ch

	output:
	set val( lg ), val( ccfpair ), file( "*ccf_prep.tsv.gz" ) into ( ccf_prep_ch )

	script:
	"""
	vcftools \
		--gzvcf ${vcf} \
		--chr ${lg} \
		--positions ${pos} \
		--recode \
		--stdout | \
		gzip > input.vcf.gz

	Rscript --vanilla \$BASE_DIR/R/ccf_prep.R \$BASE_DIR/2_analysis/geva/ ${lg} ${ccfpair.target} ${ccfpair.querry}
	"""
}

// git 12.10
// run ccf
process run_ccf {
	label 'L_20g6h_run_ccf'
	publishDir "../../2_analysis/ccf/", mode: 'copy'

	input:
	//set val( ccfidx ), val( lg ), file( ccf ) from ( ccf_prep_ch )
	set val( lg ), val( ccfpair ), file( ccf ) from ( ccf_prep_ch )

	output:
	 file( "ccf.*.gz" ) into ( ccf_output )

	script:
	"""
  # x=\$((${ccfidx}+2))

	for idx in \$( seq 1 6 ); do
		x=\$((\${idx}+2))
	  zcat ${ccf} | \
		   cut -f \$x | \
			  tail -n 2+ > input.\$x.txt

	echo "idx_"\$x > output.\$x.tsv
	ccf input.\$x.txt >> output.\$x.tsv
	done

	zcat ${ccf} | \
		cut -f -2,9 | \
		paste -d "\\t" - output.* | \
		gzip > ccf.${lg}.${ccfpair.target}.${ccfpair.querry}.tsv.gz

	"""
}
