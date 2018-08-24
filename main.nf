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
    input:
    val x from samples_ch

    output:
    file "sample_line_${x.id}.txt" into sample_lines

    script:
    """
    echo "${x.label}\t${x.file_fwd}\t${x.file_rev}\t${x.flowcell_id_fwd}\t${x.lane_fwd}\t${x.flowcell_id_fwd}.${x.lane_fwd}\t${x.company}" > sample_line_${x.id}.txt
    """
}

process bar {
  echo true
  input:
  file "sample_line_*.txt" from sample_lines.collect()

  """
  echo -e "Label\tFwd\tRev\tFlowcell\tLane\tRead_group\tCompany" > out.csv
  cat *.txt >> out.csv

  head out.csv
  """
}