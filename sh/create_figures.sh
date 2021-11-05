#/usr/bin/bash
# git 20

# preparation for Fig. 1
Rscript --vanilla R/fst_permutation.R \
    2_analysis/fst_signif/random/

# Main Figures
# -------------------------------------------------
# git 20.1
# needs to be run interactively because fonts are not accessible for Rscript
# Rscript --vanilla R/fig/plot_F1.R \
#  2_analysis/fst/50k/ \
#  2_analysis/summaries/fst_globals.txt \
#  2_analysis/summaries/fst_permutation_summary.tsv \
#  2_analysis/fotl/concat_R24ed.treefile

# git 20.2
Rscript --vanilla R/fig/plot_F2.R \
    2_analysis/fst/50k/multi_fst.50k.tsv.gz \
    2_analysis/GxP/50000/ \
    2_analysis/summaries/fst_outliers_998.tsv \
    https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
    2_analysis/twisst/weights/ \
    ressources/plugin/trees/ \
    2_analysis/summaries/fst_globals.txt

# git 20.3
Rscript --vanilla R/fig/plot_F3.R \
    2_analysis/msmc/output/ \
    2_analysis/cross_coalescence/output/ \
    2_analysis/msmc/setup/msmc_grouping.txt \
    2_analysis/msmc/setup/msmc_cc_grouping.txt \
    2_analysis/summaries/fst_globals.txt

# git 20.4
Rscript --vanilla R/fig/plot_F4.R \
    2_analysis/astral/astral_5000x_5kb_v1_noS.tre \
    2_analysis/ibd/cM_converted/no_outgr_bed95_8.conv_filterd.tsv

# git 20.5
Rscript --vanilla R/fig/plot_F5.R \
    2_analysis/twisst/weights/ \
    ressources/plugin/trees/ \
    https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
    2_analysis/summaries/fst_outliers_998.tsv \
    2_analysis/dxy/50k/ \
    2_analysis/fst/50k/ \
    2_analysis/summaries/fst_globals.txt \
    2_analysis/GxP/50000/ \
    200 \
    5 \
    2_analysis/revPoMo/outlier_regions/

# git 20.6
Rscript --vanilla R/fig/plot_F6.R \
    2_analysis/summaries/fst_outliers_998.tsv \
    2_analysis/geva/ \
    2_analysis/GxP/bySNP/

# Suppl. Figures
# -------------------------------------------------
# git 20.7
Rscript --vanilla R/fig/plot_SF1.R \
    ressources/Rabosky_etal_2018/dataFiles/ratemat_enhanced.csv

# git 20.8
# needs to be run interactively because fonts are not accessible for Rscript
# Rscript --vanilla R/fig/plot_SF2.R \
#    ressources/Rabosky_etal_2018/

# git 20.9
Rscript --vanilla R/fig/plot_SF3.R \
    2_analysis/pca/

# git 20.10
Rscript --vanilla R/fig/plot_SF4.R \
    2_analysis/fst/50k/ \
    2_analysis/summaries/fst_outliers_998.tsv \
    2_analysis/summaries/fst_globals.txt

# git 20.11
Rscript --vanilla R/fig/plot_SF5.R \
    2_analysis/fst/50k/ \
    2_analysis/summaries/fst_globals.txt

# git 20.12
Rscript --vanilla R/fig/plot_SF6.R \
    2_analysis/summaries/fst_globals.txt \
    2_analysis/fst/50k/ \
    2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz

# git 20.13
Rscript --vanilla R/fig/plot_SF7.R \
    2_analysis/dxy/50k/

# git 20.14
Rscript --vanilla R/fig/plot_SF8.R \
    2_analysis/dxy/50k/ \
    2_analysis/fst/50k/multi_fst.50k.tsv.gz \
    2_analysis/GxP/50000/ \
    2_analysis/summaries/fst_outliers_998.tsv \
    https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
    2_analysis/twisst/weights/ \
    ressources/plugin/trees/ \
    2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz \
    2_analysis/summaries/fst_globals.txt

# git 20.15
Rscript --vanilla R/fig/plot_SF9.R \
    2_analysis/pi/50k/ \
    2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz

# git 20.16
Rscript --vanilla R/fig/plot_SF10.R \
    2_analysis/dxy/50k/

# git 20.17
Rscript --vanilla R/fig/plot_SF11.R \
    2_analysis/newhyb/nh_input/NH.Results/

# git 20.18
Rscript --vanilla R/fig/plot_SF12.R \
    ressources/species_order_alpha.txt \
    2_analysis/dstats/hyp_ld05_dtrios_BBAA.txt \
    2_analysis/dstats/BBAA_ld05.csv \
    2_analysis/dstats/BBAA_sign_ld05.csv

# git 20.19
Rscript --vanilla R/fig/plot_SF13.R \
    2_analysis/summaries/fst_outliers_998.tsv

# git 20.20
Rscript --vanilla R/fig/plot_SF14.R \
    2_analysis/raxml/lg04.1_155N.raxml.support \
    2_analysis/raxml/lg12.3_155N.raxml.support \
    2_analysis/raxml/lg12.4_155N.raxml.support

# git 20.21
Rscript --vanilla R/fig/plot_SF15.R \
    2_analysis/raxml/lg04.1_hySN.raxml.support \
    2_analysis/raxml/lg12.3_hySN.raxml.support \
    2_analysis/raxml/lg12.4_hySN.raxml.support

# git 20.22
Rscript --vanilla R/fig/plot_SF16.R \
    2_analysis/admixture/ \
    metadata/phenotypes.sc

# git 20.23
Rscript --vanilla R/fig/plot_SF17.R \
    2_analysis/astral/astral_5000x_5kb_v1_all.tre

# git 20.24
Rscript --vanilla R/fig/plot_SF18.R \
    2_analysis/summaries/fst_outliers_998.tsv

# git 20.25
Rscript --vanilla R/fig/plot_SF19.R \
    2_analysis/summaries/fst_outliers_998.tsv \
    2_analysis/geva/ \
    2_analysis/GxP/bySNP/

# git 20.26
Rscript --vanilla R/fig/plot_SF20.R \
    2_analysis/GxP/50000/

# git 20.27
Rscript --vanilla R/fig/plot_SF21.R \
    2_analysis/pi/50k/

# -------------------------------------------------
rm Rplots.pdf