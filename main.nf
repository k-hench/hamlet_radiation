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

     gatk --java-options "-Xmx20G" \\
         FastqToSam \\
         SAMPLE_NAME=${x.label} \\
     		 FASTQ=\$BASE_DIR/${x.file_fwd} \\
         FASTQ2=\$BASE_DIR/${x.file_rev} \\
     		 OUTPUT=${x.label}.${x.lane_fwd}.ubam.bam \\
     		 READ_GROUP_NAME=${x.label}.${x.lane_fwd} \\
     		 LIBRARY_NAME=${x.label}".lib1" \\
     		 PLATFORM_UNIT=${x.flowcell_id_fwd}.${x.lane_fwd} \\
     		 PLATFORM=Illumina \\
     		 SEQUENCING_CENTER=${x.company} \\
     		 TMP_DIR=\$BASE_DIR/temp_files;
     """
 }