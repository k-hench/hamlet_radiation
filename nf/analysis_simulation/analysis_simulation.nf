#!/usr/bin/env nextflow
/* nextflow run 2.4.1.msms.nf -with-dag msms.png -c nextflow.config -resume */

extend = Channel.from( 500000 )
recom = Channel.from( 0.02 )
sel  = Channel.from( 0.5 )
ne = Channel.from( 1000, 10000, 100000, 100000 )
mig  = Channel.from( 0.01, 0.001, 0.0001, 0.00001)
divt  = Channel.from( 1000, 10000, 100000, 1000000 )
dom = Channel.from( 0.5 )
mode = Channel.from( 'seq', 'poly' )

par = extend
	.combine( recom )
	.combine( ne )
	.combine( sel )
	.combine( mig )
	.combine( divt )
	.combine( dom )
	.combine( mode )
	.map{ row -> [ ext:row[0], rec:row[1], ne:row[2], sel:row[3], mig:row[4], divt:row[5], dom:row[6], mode:row[7]] }

process msms_seq {
	label "L_20g2h_msms"
	tag "msms_${x.mode}_Ne_${x.ne}_mig_${x.mig}_t_${x.divt}"

	input:
	val( x ) from par

	output:
	set val( x ), file( "*.vcf.gz" ), file( "*_${x.dom}.txt" ) into ( msms_output, vcf2pca )
	set val( x ), file( "fst_decay.msms.seq_gen.fa" ) into ( vcf2geno )

	script:
	"""
	module load openssl1.0.2

	rec=\$(awk "BEGIN {print 4*${x.ne}*${x.rec}}")
	mig=\$(awk "BEGIN {print 4*${x.ne}*${x.mig}}")


	# check simulation scenario
	if [ "${x.mode}" == "seq" ];then

		sel1=\$(awk "BEGIN {print 2*${x.ne}*${x.sel}}") # sequential
		sel2=\$(awk "BEGIN {print 2*${x.ne}*${x.sel}}")
		sel3=\$(awk "BEGIN {print 2*${x.ne}*${x.sel}}")
		sel4=\$(awk "BEGIN {print 2*${x.ne}*${x.sel}}")

		divt1=\$(awk "BEGIN {print ${x.divt}/(4*${x.ne}*1)}") # sequential
		divt2=\$(awk "BEGIN {print ${x.divt}/(4*${x.ne}*2)}")
		divt3=\$(awk "BEGIN {print ${x.divt}/(4*${x.ne}*4)}")

	else

		sel1=\$(awk "BEGIN {print 2*${x.ne}*${x.sel}/1}") # politomy
		sel2=\$(awk "BEGIN {print 2*${x.ne}*${x.sel}/2}")
		sel3=\$(awk "BEGIN {print 2*${x.ne}*${x.sel}*0}")
		sel4=\$(awk "BEGIN {print 2*${x.ne}*${x.sel}/1}")

		divt1=\$(awk "BEGIN {print ${x.divt}/(4*${x.ne})}") # politomy
		divt2=\$(awk "BEGIN {print ${x.divt}/(4*${x.ne})}")
		divt3=\$(awk "BEGIN {print ${x.divt}/(4*${x.ne})}")

	fi

	dom1=\$(awk -v s=\$sel1 "BEGIN {print s*${x.dom}}")
	dom2=\$(awk -v s=\$sel2 "BEGIN {print s*${x.dom}}")
	dom3=\$(awk -v s=\$sel3 "BEGIN {print s*${x.dom}}")
	dom4=\$(awk -v s=\$sel4 "BEGIN {print s*${x.dom}}")
	antidom1=\$(awk -v s=\$sel1 "BEGIN {print s*(1-${x.dom})}")
	antidom2=\$(awk -v s=\$sel2 "BEGIN {print s*(1-${x.dom})}")
	antidom3=\$(awk -v s=\$sel3 "BEGIN {print s*(1-${x.dom})}")
	antidom4=\$(awk -v s=\$sel4 "BEGIN {print s*(1-${x.dom})}")

	msms 48 1 \
		-T \
		-I 4 12 12 12 12 \
		-ej \$divt1 4 1 \
		-ej \$divt2 3 1 \
		-ej \$divt3 2 1 \
		-m 1 4 \$mig \
		-m 4 1 \$mig \
		-m 2 4 \$mig \
		-m 4 2 \$mig \
		-m 3 4 \$mig \
		-m 4 3 \$mig \
		-m 3 2 \$mig \
		-m 2 3 \$mig \
		-m 3 1 \$mig \
		-m 1 3 \$mig \
		-m 2 1 \$mig \
		-m 1 2 \$mig \
		-r \$rec ${x.ext} \
		-Sc 0 1 0 \$antidom4 \$sel4 \
		-Sc 0 2 \$sel2 \$dom2 0 \
		-Sc 0 3 \$sel3 \$dom3 0 \
		-Sc 0 4 \$sel1 \$dom1 0 \
		-Sp 0.5 \
		-SI \$divt1 4 0.5 0.5 0.5 0.5 \
		-N ${x.ne} \
		-threads 4 \
		-seed 27678 | \
		grep ";"  > fst_decay.msms.grep.${x.ext}_${x.rec}_${x.ne}_${x.sel}_${x.mig}_${x.divt}_${x.dom}.txt


	partitions=\$(wc -l fst_decay.msms.grep.${x.ext}_${x.rec}_${x.ne}_${x.sel}_${x.mig}_${x.divt}_${x.dom}.txt)

	seq-gen -of -mHKY -l ${x.ext} -s 0.01 -p \$partitions < fst_decay.msms.grep.${x.ext}_${x.rec}_${x.ne}_${x.sel}_${x.mig}_${x.divt}_${x.dom}.txt  | \
		sed 's/^>\\([0-9]\\)\$/>0\\1/g' > fst_decay.msms.seq_gen.fa

	cat fst_decay.msms.seq_gen.fa | msa2vcf > fst_decay.seq-gen.${x.ext}_${x.rec}_${x.ne}_${x.sel}_${x.mig}_${x.divt}_${x.dom}.vcf
	bgzip fst_decay.seq-gen.${x.ext}_${x.rec}_${x.ne}_${x.sel}_${x.mig}_${x.divt}_${x.dom}.vcf
	"""
}

process prep_dxy {
	label "L_20g2h_msms"
	tag "prep_${x.mode}_Ne_${x.ne}_mig_${x.mig}_t_${x.divt}"

	input:
	set val( x ), file( fa ) from vcf2geno

	output:
	set val( x ), file( "run.${x.ext}_${x.rec}_${x.ne}_${x.sel}_${x.mig}_${x.divt}_${x.dom}.geno.gz" ) into ( dxy_prep )
	set val( x ), file( "vcf.gz" ) into msms_vcf

	script:
	"""
	module load openssl1.0.2

	cat ${fa} | msa2vcf -a  | bgzip > vcf.gz

	python \$SFTWR/genomics_general/VCF_processing/parseVCF.py -i vcf.gz |  \
	bgzip > run.${x.ext}_${x.rec}_${x.ne}_${x.sel}_${x.mig}_${x.divt}_${x.dom}.geno.gz
	"""
}

process dxy_run {
	label "L_20g2h_msms"
	tag "dxy_${x.mode}_Ne_${x.ne}_mig_${x.mig}_t_${x.divt}"
	publishDir "../../2_analysis/simulation/dxy", mode: 'copy'

	input:
	set val( x ), file( geno ) from ( dxy_prep )

	output:
	file( "run.*.dxy.gz" ) into ( dxy_output )

	script:
	"""
	for k in {01..12};do echo -e "\$k\tA" >> pop.txt; done
	for k in {13..24};do echo -e "\$k\tB" >> pop.txt; done
	for k in {25..36};do echo -e "\$k\tC" >> pop.txt; done
	for k in {37..48};do echo -e "\$k\tD" >> pop.txt; done

	for k in A-B A-C A-D B-C B-D C-D; do
		pop1=\$(echo \$k | cut -c 1)
		pop2=\$(echo \$k | cut -c 3)

		python \$SFTWR/genomics_general/popgenWindows.py \
			-w 10000 -s 1000 \
			--popsFile pop.txt \
			-p \$pop1 -p \$pop2 \
			-g ${geno} \
			-o run.${x.ext}_${x.rec}_${x.ne}_${x.sel}_${x.mig}_${x.divt}_${x.dom}_\$pop1-\$pop2-${x.mode}.dxy.gz \
			-f phased --writeFailedWindows -T 1
	done
	"""
}

process fst_run {
	label "L_20g2h_msms"
	tag "fst_${x.mode}_Ne_${x.ne}_mig_${x.mig}_t_${x.divt}"
	publishDir "../../2_analysis/simulation/fst", mode: 'copy'

	input:
	set val( x ), file( vcf ) from  msms_vcf

	output:
	file( "run.*.windowed.weir.fst" ) into fst_50k_output

	script:
	"""
	for k in {01..12};do echo \$k >> popA.txt; done
	for k in {13..24};do echo \$k >> popB.txt; done
	for k in {25..36};do echo \$k >> popC.txt; done
	for k in {37..48};do echo \$k >> popD.txt; done

	for k in A-B A-C A-D B-C B-D C-D; do

		pop1=\$(echo \$k | cut -c 1)
		pop2=\$(echo \$k | cut -c 3)

		vcftools --gzvcf ${vcf} \
			--min-alleles 2 \
			--weir-fst-pop pop\$pop1.txt \
			--weir-fst-pop pop\$pop2.txt \
			--fst-window-step 1000 \
			--fst-window-size 10000 \
			--out run.${x.ext}_${x.rec}_${x.ne}_${x.sel}_${x.mig}_${x.divt}_${x.dom}_\$pop1-\$pop2-${x.mode}
	done
	"""
}
