alias nf_run_gatk="nextflow run genotyping.nf -with-dag docs/genotyping.dot -c nextflow.config -resume"
alias nf_run_basic="nextflow run analysis_pca_admx_fst_gxp.nf -with-dag docs/analysis_basic.dot -c nextflow.config -resume"
alias nf_run_phylo="nextflow run analysis_fasttree_twisst.nf -with-dag docs/analysis_phylo.dot -c nextflow.config -resume"
alias nf_run_msmc="nextflow run analysis_msmc.nf -with-dag docs/analysis_msmc.dot -c nextflow.config -resume"