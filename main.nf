#!/usr/bin/env nextflow

 /* open the pipeline based on the metadata spread sheet that includes all
  * information necessary to assign read groups to the sequencing data */
 params.index = 'metadata/file_info.txt'

 /* split the spread sheet by row and feed it into a channel */
 Channel
     .fromPath(params.index)
     .splitCsv(header:true, sep:"\t")
     .map{ row -> [ id:row.id, label:row.label, file_fwd:row.file_fwd, file_rev:row.file_rev, flowcell_id_fwd:row.flowcell_id_fwd, lane_fwd:row.lane_fwd, company:row.company] }
     .set { samples_ch }

 /* for every sequencing file, convert into ubam format and assign read groups */
 process split_samples {
     label 'L_20g2h_split_samples'
     conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'

     input:
     val x from samples_ch

     output:
     set val( "${x.label}.${x.lane_fwd}" ), file( "${x.label}.${x.lane_fwd}.ubam.bam" ) into ubams_mark, ubams_merge

     script:
     """
     echo -e "---------------------------------"
     echo -e "Label:\t\t${x.label}\nFwd:\t\t${x.file_fwd}\nRev:\t\t${x.file_rev}"
     echo -e "Flowcell:\t${x.flowcell_id_fwd}\nLane:\t\t${x.lane_fwd}"
     echo -e "Read group:\t${x.flowcell_id_fwd}.${x.lane_fwd}\nCompany:\t${x.company}"

     mkdir -p \$BASE_DIR/temp_files

     gatk --java-options "-Xmx20G" \
         FastqToSam \
         -SM=${x.label} \
     		 -F1=\$BASE_DIR/data/seqdata/${x.file_fwd} \
         -F2=\$BASE_DIR/data/seqdata/${x.file_rev} \
     		 -O=${x.label}.${x.lane_fwd}.ubam.bam \
     		 -RG=${x.label}.${x.lane_fwd} \
     		 -LB=${x.label}".lib1" \
     		 -PU=${x.flowcell_id_fwd}.${x.lane_fwd} \
     		 -PL=Illumina \
     		 -CN=${x.company} \
     		 --TMP_DIR=\$BASE_DIR/temp_files;
     """
 }

 /* for every ubam file, mark Illumina adapters */
 process mark_adapters {
   label 'L_20g2h_mark_adapters'
   conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
   tag "${sample}"

   input:
   set val( sample ), file( input ) from ubams_mark

   output:
   set val( sample ), file( "*.adapter.bam") into adapter_bams
   file "*.adapter.metrics.txt" into adapter_metrics

   script:
   """
  gatk --java-options "-Xmx18G" \
        MarkIlluminaAdapters \
        -I=${input} \
        -O=${sample}.adapter.bam \
        -M=${sample}.adapter.metrics.txt \
        -TMP_DIR=\$BASE_DIR/temp_files;
   """
 }

 adapter_bams
     .combine(ubams_merge, by:0)
     .set {merge_input}

 /* this step includes a 3 step pipeline:
  *  - re-transformatikon into fq format
  *  - mapping aginst the reference genome_file
  *  - merging with the basuch ubams to include
       read group information */

 process map_and_merge {
   label 'L_75g24h8t_map_and_merge'
   conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
   tag "${sample}"

   input:
   set val( sample ), file( adapter_bam_input ), file( ubam_input ) from merge_input

   output:
   set val( sample ), file( "*.mapped.bam" ) into mapped_bams

   script:
   """
   set -o pipefail
   gatk --java-options "-Xmx68G" \
        SamToFastq \
        -I=${adapter_bam_input} \
        -FASTQ=/dev/stdout \
        -INTERLEAVE=true \
        -NON_PF=true \
        -TMP_DIR=\$BASE_DIR/temp_files | \
    bwa mem -M -t 8 -p \$BASE_DIR/ressources/HP_genome_unmasked_01.fa /dev/stdin |
    gatk --java-options "-Xmx68G" \
        MergeBamAlignment \
         --VALIDATION_STRINGENCY SILENT \
         --EXPECTED_ORIENTATIONS FR \
         --ATTRIBUTES_TO_RETAIN X0 \
         -ALIGNED_BAM=/dev/stdin \
         -UNMAPPED_BAM=${ubam_input} \
         -OUTPUT=${sample}.mapped.bam \
         --REFERENCE_SEQUENCE=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz \
         -PAIRED_RUN true \
         --SORT_ORDER "unsorted" \
         --IS_BISULFITE_SEQUENCE false \
         --ALIGNED_READS_ONLY false \
         --CLIP_ADAPTERS false \
         --MAX_RECORDS_IN_RAM 2000000 \
         --ADD_MATE_CIGAR true \
         --MAX_INSERTIONS_OR_DELETIONS -1 \
         --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
         --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
         --ALIGNER_PROPER_PAIR_FLAGS true \
         --UNMAP_CONTAMINANT_READS true \
         -TMP_DIR=\$BASE_DIR/temp_files
   """
 }

 /* for every mapped sample,sort and mark duplicates
 * (intermediate step is required to create .bai file) */
 process mark_duplicates {
   label 'L_32g30h_mark_duplicates'
   conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
   publishDir "1_genotyping/0_sorted_bams/", mode: 'symlink'
   tag "${sample}"

   input:
   set val( sample ), file( input ) from mapped_bams

   output:
   set val { sample  - ~/\.(\d+)/ }, val( sample ), file( "*.dedup.bam") into dedup_bams
   file "*.dedup.metrics.txt" into dedup_metrics

   script:
   """
   set -o pipefail
   gatk --java-options "-Xmx30G" \
        SortSam \
        -I=${input} \
        -O=/dev/stdout \
        --SORT_ORDER="coordinate" \
        --CREATE_INDEX=false \
        --CREATE_MD5_FILE=false \
        -TMP_DIR=\$BASE_DIR/temp_files \
        | \
  gatk --java-options "-Xmx30G" \
      SetNmAndUqTags \
      --INPUT=/dev/stdin \
      --OUTPUT=intermediate.bam \
      --CREATE_INDEX=true \
      --CREATE_MD5_FILE=true \
      -TMP_DIR=\$BASE_DIR/temp_files \
      --REFERENCE_SEQUENCE=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz

   gatk --java-options "-Xmx30G" \
        MarkDuplicates \
        -I=intermediate.bam \
        -O=${sample}.dedup.bam \
        -M=${sample}.dedup.metrics.txt \
        -MAX_FILE_HANDLES=1000  \
        -TMP_DIR=\$BASE_DIR/temp_files

   rm intermediate*
   """
 }

 /* index al bam files */
process index_bam {
  label 'L_32g1h_index_bam'
  conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
  tag "${sample}"

  input:
  set val( sample ), val( sample_lane ), file( input ) from dedup_bams

  output:
  set val( sample ), val( sample_lane ), file( input ), file( "*.bai") into indexed_bams

  script:
  """
  gatk --java-options "-Xmx30G" \
  BuildBamIndex \
  -INPUT=${input}
  """
}

 /* collect all bam files for each sample */
indexed_bams
  .groupTuple()
  .set {tubbled}

 /* create one *.g.vcf file per sample */
process receive_tuple {
  label 'L_36g47h_receive_tuple'
  conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
  tag "${sample}"

  input:
  set sample, sample_lane, bam, bai from tubbled

  output:
  file( "*.g.vcf.gz") into gvcfs
  file( "*.vcf.gz.tbi") into tbis

  script:
  """
  INPUT=\$(echo ${bam}  | sed  's/\\[/-I /g; s/\\]//g; s/,/ -I/g')

  gatk --java-options "-Xmx35g" HaplotypeCaller  \
    -R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
    \$INPUT \
    -O ${sample}.g.vcf.gz \
    -ERC GVCF
  """
}

 /* collect and combine all *.g.vcf files */
process gather_gvcfs {
  label 'L_O88g90h_gather_gvcfs'
  conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
  echo true

  input:
  file( gvcf ) from gvcfs.collect()
  file( tbi ) from tbis.collect()

  output:
  set file( "cohort.g.vcf.gz" ), file( "cohort.g.vcf.gz.tbi" ) into gcvf_snps, gvcf_acs

  script:
  """
  GVCF=\$(echo " ${gvcf}" | sed 's/ /-V /g; s/vcf.gz/vcf.gz /g')

  gatk --java-options "-Xmx85g" \
  CombineGVCFs \
  -R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
  \$GVCF \
  -O cohort.g.vcf.gz
  """
}

/* actual genotyping step (varinat sites only) */
process joint_genotype_snps {
  label 'L_O88g90h_joint_genotype'
  conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
  publishDir "1_genotyping/1_raw_vcfs/", mode: 'symlink'

  input:
  set file( vcf ), file( tbi ) from gcvf_snps

  output:
  set file( "raw_var_sites.vcf.gz" ), file( "raw_var_sites.vcf.gz.tbi" ) into ( raw_var_sites, raw_var_sites_to_metrics )

  script:
  """
  gatk --java-options "-Xmx85g" \
  GenotypeGVCFs \
  -R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
  -V=${vcf} \
  -O=intermediate.vcf.gz

  gatk --java-options "-Xmx85G" \
  SelectVariants \
  -R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
  -V=intermediate.vcf.gz \
  --select-type-to-include=SNP \
  -O=raw_var_sites.vcf.gz

  rm intermediate.*
  """
}

/* generate a LG channel */
/*
LG_ids1 = Channel.from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
*/
/* actual genotyping step
 * (all callable sites,
 *  one process per LG) */
/*
process joint_genotype_acs {
  label 'L_105g30h_joint_genotype_acs'
  conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
  publishDir "1_genotyping/1_raw_vcfs/", mode: 'symlink'

  input:
  set file( vcf ), file( tbi ) from gvcf_acs
  val LGid from LG_ids1

  output:
  file( "raw_gvcf_acs.*" ) into raw_acs_by_ls

  script:
  """
  gatk --java-options "-Xmx100g" \
  GenotypeGVCFs \
  --includeNonVariantSites \
  -R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
  -V ${vcf} \
  -L LG${LGid} \
  -O intermediate.vcf.gz

  gatk --java-options "-Xmx100G" \
  SelectVariants \
  -R \$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
  -V intermediate.vcf.gz \
  --select-type-to-include=SNP \
  -O=raw_gvcf_acs.LG${LGid}.vcf.gz

  rm intermediate.*
  """
}
*/
/* produce metrics table to determine filtering thresholds - ups forgot to extract SNPS first*/
process joint_genotype_metrics {
  label 'L_28g5h_genotype_metrics'
  conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
  publishDir "1_genotyping/1_raw_vcfs/", mode: 'move'

  input:
  set file( vcf ), file( tbi ) from raw_var_sites_to_metrics

  output:
  file( "${vcf}.table.txt" ) into raw_metrics

  script:
  """
  gatk --java-options "-Xmx25G" \
  VariantsToTable \
  --variant=${vcf} \
  --output=${vcf}.table.txt \
  -F=CHROM -F=POS -F=MQ \
  -F=QD -F=FS -F=MQRankSum -F=ReadPosRankSum \
  --show-filtered
  """
}


/* produce metrics table to determine filtering thresholds */
process filterSNPs {
  label 'L_78g10h_filter_Snps'
  publishDir "1_genotyping/2_gatk_filtered/", mode: 'symlink'

  input:
  set file( vcf ), file( tbi ) from raw_var_sites

  output:
  set file( "filterd_bi-allelic.vcf.gz" ), file( "filterd_bi-allelic.vcf.gz.tbi" ) into filtered_snps

  script:
  """
  gatk --java-options "-Xmx75G" \
		VariantFiltration \
		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
		-V ${vcf} \
		-O=intermediate.vcf.gz \
		--filter-expression "QD < 2.5" \
		--filter-name "filter_QD" \
		--filter-expression "FS > 25.0" \
		--filter-name "filter_FS" \
		--filter-expression "MQ < 52 || MQ > 65" \
		--filter-name "filter_MQ" \
		--filter-expression "MQRankSum < -0.2 || MQRankSum > 0.2" \
		--filter-name "filter_MQRankSum" \
		--filter-expression "ReadPosRankSum < -2 || ReadPosRankSum > 2 " \
		--filter-name "filter_ReadPosRankSum"

		gatk --java-options "-Xmx75G" \
	  SelectVariants \
	  -R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \
	  -V=intermediate.vcf.gz \
		-O=intermediate.filterd.vcf.gz \
		--select-type-to-include=SNP
		--exclude-filtered

		vcftools \
			--gzvcf intermediate.filterd.vcf.gz \
			--max-missing-count 17 \
			--max-alleles 2 \
			--stdout  \
			--recode | \
			bgzip > filterd_bi-allelic.vcf.gz

		tabix -p vcf filterd_bi-allelic.vcf.gz

	  rm intermediate.*
  """
}


