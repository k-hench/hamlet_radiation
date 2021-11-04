#!/usr/bin/env nextflow

// ----------------------- DISCLAIMER ----------------------
// this pipeline was not actually run using nextflow,
// but managed manually
// ---------------------------------------------------------

// git 19.1
Channel
	.fromPath('../../ressources/serraninae/Serraninae_Rabosky.phy')
	.set{ serraninae_phy_ch }

// git 19.2
Channel
	.fromPath('../../ressources/serraninae/partitions.txt')
	.set{ partitions_ch }

// git 19.3
// Obtain gene boundaries and split alignment into individual genes
process get_gene_boundaries  {
	label 'gene_boundaries'

	input:
	set file( phy ), file( part ) from serraninae_phy_ch.combine( partitions_ch )

	output:
	file( "*serrR.aln" ) into ( gene_boundaries_ch, gene_boundaries_ch2 )

	script:
	"""
	# manually extract Serraninae from Rabosky alignment, final_alignment.phylip (FToL Dryad) > Serraninae_Rabosky.phy
	perl \$BASE_DIR/pl/Phylip2Fasta.pl ${phy} Serraninae_Rabosky.aln

	sed 's/-/N/g' Serraninae_Rabosky.aln > Serraninae_RaboskyN.aln

	sed 's/DNA, \\(.*\\) = \\(.*\\)-\\(.*\\)/\\1\\t\\2\\t\\3/g' ${part} > partitions.bed   # (partitions.txt from FToL Dryad)

	awk '{ print \$2 }' partitions.bed | tail -n +2 | tr '\\n' ',' | sed 's/.\$//'

	# phast from http://compgen.cshl.edu/phast/
	\$SFTWR/phast/bin/msa_split Serraninae_RaboskyN.aln -r gene \
	  --by-index 980,1757,2292,2974,4115,4955,5681,6614,7244,7988,8822,9869,11402,12142,12952,13663,15117,16341,17265,17910,18633,19728,20715,21522,22350,23115

	awk '{ print \$1, \$2, \$3 }' partitions.bed | xargs -n 3 sh -c 'mv gene.\$1-\$2.fa \$0_serrR.aln'
	"""
}

// git 19.4
// Identify genes in reference genome
process identify_genes  {
	label 'identify_genes'

	input:
	file( aln ) from gene_boundaries_ch

	output:
	file( "coord_R24_Hpue_ed.bed" ) into gene_coords_ch

	script:
	"""
	# manually select query sequences from Serraninae_RaboskyN.aln > queries_R24.fas

	mkdir blast_db
	cp \$BASE_DIR/ressources/HP_genome_unmasked_01.fa blast_db/Hpue_genome_unmasked_01.fas
	makeblastdb -in blast_db/Hpue_genome_unmasked_01.fas -dbtype nucl -parse_seqids

	blastn -query queries_R24.fas -db blast_db/Hpue_genome_unmasked_01.fas -out blast_R24-Hpue_aln.txt -outfmt 0 -evalue 1e-10
	blastn -query queries_R24.fas -db blast_db/Hpue_genome_unmasked_01.fas -out Rblast_R24-Hpue_tab.csv -outfmt 6 -evalue 1e-10

	awk 'OFS = "\\t" { if (\$9 < \$10) print \$2, \$9, \$10, \$1, "+" ; else print \$2, \$10, \$9, \$1, "-" }' blast_R24-Hpue_tab.csv | \
		sed 's/\\(.*\\)_\\(.*\\)\\t/\\1\\t/g' > coord_R24_Hpue.csv

	# manually combine segmented genes (16s, glyt, tbr1, zic1) in coord_R24_Hpue.csv > coord_R24_Hpue_ed.csv

	sed -Ei 's/_[A-Z][a-z]{3}//g' coord_R24_Hpue_ed.bed
	"""
}

// git 19.5
Channel
	.fromPath('../../1_genotyping/2_raw_vcfs/all_sites.unplaced.vcf.gz')
	.set{ vcf_unlplaced_ch }

// git 19.6
// Prepare contigs 
process prepare_contigs  {
	label 'prepare_contigs'

	input:
	file ( vcf ) from vcf_unlplaced_ch

	output:
	set file( "all_sites.Contig*.vcf.gz" ), file( "*.vcf.gz.tbi" ) into contigs_ch

	script:
	"""
	zcat < ${vcf} | grep -e '#' -e 'Contig11544' | bgzip > all_sites.Contig11544.vcf.gz
	zcat < ${vcf} | grep -e '#' -e 'Contig11607' | bgzip > all_sites.Contig11607.vcf.gz
	zcat < ${vcf} | grep -e '#' -e 'Contig11888' | bgzip > all_sites.Contig11888.vcf.gz

	tabix -p vcf *.vcf.gz
	"""
}

// git 19.7
Channel
	.fromPath('../../ressources/serraninae/samples_to_include.ids')
	.set{ samples_ch }

// git 19.8
// Extract genotypes and convert to Fasta
process extract_genotypes  {
	label 'extract_genotypes'

	input:
	set file( crds ), file ( samples ) from gene_coords_ch.combine( samples_ch )

	output:
	file( "*.fas" ) into genotypes_ch

	script:
	"""
	# select representative samples based on highest coverage (except outliers) and regional diversity > samples_to_include.ids

	awk '{ print $1, $2, $3, $4 }' ${crds} | \
		xargs -n 4 sh -c 'vcftools --gzvcf \$BASE_DIR/2_raw_vcfs/all_sites."\$0".vcf.gz --keep ${samples} --chr "\$0" --from-bp "\$1" --to-bp "\$2" --recode --stdout | grep -v "##" > "\$3"_hypS.vcf'

	for FILE in ./*.vcf
		do
		vcf-to-tab < \$FILE > \${FILE%.vcf}.tab
		perl \$SFTWR/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i \${FILE%.vcf}.tab > \${FILE%.vcf}.fas
		rm \${FILE%.vcf}.tab*
	done

	rev=`awk '\$5 == "-" { print \$4 }' ${crds}`

	for GENE in \$rev
		do
		seqtk seq -r \${GENE}_hypS.fas > \${GENE}_hypS_rev.fas
		mv \${GENE}_hypS.fas \${GENE}_hypS.fas_org
	done
	"""
}

// git 19.9
Channel
	.fromPath('../../ressources/serraninae/taxa_to_exclude.ids')
	.set{ taxa_ch }

// git 19.10
// Combine Rabosky and present data
process combine_data  {
	label 'combine_data'

	input:
	set file( aln ), file( taxa ), file( fa ) from gene_boundaries_ch2.combine( taxa_ch ).combine( genotypes_ch )

	output:
	file( "*_new.fas" ) into data_combined_ch

	script:
	"""
	# manually prepared list of taxa to exclude from Rabosky's alignments > taxa_to_exclude.ids

	for FILE in *.aln
		do
		sed -e 's/N//g' -e 's/> />/g' -e '/^\$/d' \$FILE | 
		awk 'BEGIN { while ((getline < "${taxa}") >0) l[">"\$1]=1 } /^>/ { f=!l[\$1] }f' > \${FILE%.aln}_pre.fa
	done

	for FILE in *.fas
		do
		NAME=\${FILE%_hypS.fas}
		cat \$FILE \${NAME}_serrR_pre.fa > \${NAME}_new.fas
	done
	"""
}

// git 19.11
// Align sequences and remove ambiguously aligned positions
process align_sequences  {
	label 'align_sequences'

	input:
	file( fa ) from data_combined_ch

	output:
	file( "*.aln" ) into align_sequences_ch

	script:
	"""
	for FILE in ./*_new.fas
		do 
		mafft --maxiterate 1000 --globalpair --adjustdirectionaccurately \${FILE} > \${FILE%.fas}.aln
	done

	Gblocks 3_alignments/12s_new.aln -t=D -b2=0 -b3=4 -b4=5 -b5=a -e=.gb -v=100
	Gblocks 3_alignments/16s_new.aln -t=D -b2=0 -b3=4 -b4=5 -b5=a -e=.gb -v=100
	Gblocks 3_alignments/coi_new.aln -t=D -b2=0 -b3=4 -b4=5 -b5=a -e=.gb -v=100
	Gblocks 3_alignments/cytb_new.aln -t=D -b2=0 -b3=4 -b4=5 -b5=a -e=.gb -v=100

	# finalize with manual alignment editing step:
	# 12S:  removed 45±3 bp (whole gap)           > 12s_new_ed.aln
	# 16S:  removed positions 194, 252, 315       > 16s_new_ed.aln
	# coi:  removed entire gene in torpan, tabhon > coi_new_ed.aln
	# cytb: removed entire gene in torpan, tabhon > cytb_new_ed.aln
	"""
}

// ----------------------- DISCLAIMER ----------------------
// Between git 19.11 and git 19.12 the alignments 
// were additionally currated manually
// ---------------------------------------------------------

// git 19.12
// Phylogenetic inference (IQ-TREE)
process phylogenetic_inference  {
	label 'phylogenetic_inference'
	publishDir "../../2_analysis/fotl/", mode: 'copy' 

	input:
	file( aln ) from align_sequences_ch
	
	output:
	file( "concat_R24ed.treefile" ) into phylo_out

	script:
	"""
	mkdir alignments
	cp *.aln alignments/
	iqtree2 -p ./alignments --prefix concat_R24ed -B 1000 -T AUTO
	"""
}