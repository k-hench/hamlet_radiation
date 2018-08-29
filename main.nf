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
     file "${x.label}.${x.lane_fwd}.ubam.bam" into ubams

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

   input:
   file input from ubams

   output:
   file "*.adapter.bam" into adapter_bams
   file "*.adapter.metrics.txt" into adapter_metrics

   script:
   """
   BASE_FILE=\$( echo ${input} | sed 's/.ubam.bam//g' )

   gatk --java-options "-Xmx18G" \
        MarkIlluminaAdapters \
        -I=${input} \
        -O=\$BASE_FILE.adapter.bam \
        -M=\$BASE_FILE.adapter.metrics.txt \
        -TMP_DIR=\$BASE_DIR/temp_files;
   """
 }



 /* for every ubam file, mark Illumina adapters */
 process map_and_merge {
   conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'

   input:
   file input from adapter_bams

   output:
   file "*.mapped.bam" into mapped_bams

   script:
   """
   BASE_FILE=\$( echo ${input} | sed 's/.adapter.bam//g' )

  gatk --java-options "-Xmx38G" \
        SamToFastq \
        -I=${input} \
        -FASTQ=/dev/stdout \
        -CLIPPING_ATTRIBUTE=XT \
        -CLIPPING_ACTION=2 \
        -INTERLEAVE=true \
        -NON_PF=true \
        -TMP_DIR=\$BASE_DIR/temp_files | \
    bwa mem -M -t 8 -p \$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz /dev/stdin | \
    gatk --java-options "-Xmx50G" \
        -MergeBamAlignment \
        -ALIGNED_BAM=/dev/stdin \
        -UNMAPPED_BAM=$WORK/2_output/01_Ubams_BH/${BASH_REMATCH[1]}.bam \
        -OUTPUT=\$BASE_FILE.mapped.bam \
        -R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa.gz \
        -CREATE_INDEX=true \
        -ADD_MATE_CIGAR=true \
        -CLIP_ADAPTERS=false \
        -CLIP_OVERLAPPING_READS=true \
        -INCLUDE_SECONDARY_ALIGNMENTS=true \
        -MAX_INSERTIONS_OR_DELETIONS=-1 \
        -PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        -ATTRIBUTES_TO_RETAIN=XS \
        -TMP_DIR=\$BASE_DIR/temp_files;
   """
 }