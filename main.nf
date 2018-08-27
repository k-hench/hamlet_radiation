#!/usr/bin/env nextflow

/*
 *  don't forget to comment...
 */

 params.index = 'metadata/file_info.txt'

 Channel
     .fromPath(params.index)
     .splitCsv(header:true, sep:"\t")
     .map{ row -> [ id:row.id, label:row.label, file_fwd:row.file_fwd, file_rev:row.file_rev, flowcell_id_fwd:row.flowcell_id_fwd, lane_fwd:row.lane_fwd, company:row.company] }
     .set { samples_ch }

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

 process mark_adapters {
   conda '/sfs/fs6/home-geomar/smomw287/miniconda2/envs/gatk'

   input:
   file "input.ubam.bam" from ubams

   output:
   file "output.adapter.bam" into adapter_bams
   file "output.adapter.metrics.txt" into adapter_metrics

   script:
   """
   gatk --java-options ""-Xmx18G" \
        MarkIlluminaAdapters \
        -I=input.ubam.bam \
        -O=output.adapter.bam \
        -M=output.adapter.metrics.txt \
        -TMP_DIR=\$BASE_DIR/temp_files;
   """
 }