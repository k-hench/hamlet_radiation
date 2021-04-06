#/usr/bin/bash
# git 1X

Rscript --vanilla R/fig/plot_F1.R \
   2_analysis/dxy/50k/ 2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt

Rscript --vanilla R/fig/plot_F2.R \
   2_analysis/msmc/output/ 2_analysis/cross_coalescence/output/ \
   2_analysis/msmc/setup/msmc_grouping.txt 2_analysis/msmc/setup/msmc_cc_grouping.txt \
   2_analysis/summaries/fst_globals.txt

Rscript --vanilla R/fig/plot_F3.R \
   2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt

Rscript --vanilla R/fig/plot_F4.R 2_analysis/dxy/50k/ \
   2_analysis/fst/50k/multi_fst.50k.tsv.gz 2_analysis/GxP/50000/ \
   2_analysis/summaries/fst_outliers_998.tsv \
   https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
   2_analysis/twisst/weights/ ressources/plugin/trees/ \
   2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz \
   2_analysis/summaries/fst_globals.txt

Rscript --vanilla R/fig/plot_F5.R \
  2_analysis/twisst/weights/ ressources/plugin/trees/ \
  https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
  2_analysis/summaries/fst_outliers_998.tsv 2_analysis/dxy/50k/ \
  2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt \
  2_analysis/GxP/50000/ 200 5 2_analysis/fst/poptree/summary/

Rscript --vanilla R/fig/plot_SF1.R 2_analysis/newhyb/nh_input/NH.Results/

Rscript --vanilla R/fig/plot_SF2.R 2_analysis/fst/50k/ \
  2_analysis/summaries/fst_outliers_998.tsv \
  2_analysis/summaries/fst_globals.txt

Rscript --vanilla R/fig/plot_SF3.R \
  2_analysis/summaries/fst_globals.txt \
  2_analysis/fst/50k/ \
  2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz

Rscript --vanilla R/fig/plot_SF4.R 2_analysis/dxy/50k/

Rscript --vanilla R/fig/plot_SF5.R 2_analysis/pi/50k/ \
  2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz

Rscript --vanilla R/fig/plot_SF6.R 2_analysis/fasttree/ no_outgroups.all.whg.SNP.tree

Rscript --vanilla R/fig/plot_SF7.R 2_analysis/admixture/ metadata/phenotypes.sc

Rscript --vanilla R/fig/plot_SF8.R 2_analysis/pi/50k/

Rscript --vanilla R/fig/plot_SF9.R 2_analysis/GxP/50000/
