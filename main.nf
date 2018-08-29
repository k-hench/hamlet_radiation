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
   conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
   tag "${sample}"

   input:
   set val( sample ), file( adapter_bam_input ), file( ubam_input ) from merge_input

   output:
   set val( sample ), file( "*.mapped.bam" ) into mapped_bams

   script:
   """
   set -o pipefail
   gatk --java-options "-Xmx48G" \
        SamToFastq \
        -I=${adapter_bam_input} \
        -FASTQ=/dev/stdout \
        -INTERLEAVE=true \
        -NON_PF=true \
        -TMP_DIR=\$BASE_DIR/temp_files | \
    bwa mem -M -t 8 -p \$BASE_DIR/ressources/HP_genome_unmasked_01.fa /dev/stdin |
    gatk --java-options "-Xmx48G" \
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
   conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'
   publishDir "1_genotyping/0_sorted_bams/", mode: 'symlink'
   tag "${sample}"

   input:
   set val( sample ), file( input ) from mapped_bams

   output:
   set val( sample ), file( "*.dedup.bam") into dedup_bams
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