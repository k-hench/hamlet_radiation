#/usr/bin/bash

Rscript --vanilla R/fig/plot_F1.R 2_analysis/dxy/50k/ 2_analysis/fst/50k/

Rscript --vanilla R/fig/plot_F2.R 2_analysis/dxy/50k/ \
   2_analysis/fst/50k/multi_fst.50k.tsv.gz 2_analysis/GxP/50000/ \
   2_analysis/summaries/fst_outliers_998.tsv \
   https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
   2_analysis/twisst/weights/ ressources/plugin/trees/ \
   2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz \
   2_analysis/summaries/fst_globals.txt

Rscript --vanilla R/fig/plot_F3.R \
  2_analysis/twisst/weights/ ressources/plugin/trees/ \
  https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
  2_analysis/summaries/fst_outliers_998.tsv 2_analysis/dxy/50k/ \
  2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt \
  2_analysis/GxP/50000/ 200 5

Rscript --vanilla R/fig/plot_SFx1.R \
  2_analysis/summaries/fst_globals.txt \
  2_analysis/fst/50k/ \
  2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz

Rscript --vanilla R/fig/plot_SFx2.R \
  metadata/phenotypes.sc ressources/img/

Rscript --vanilla R/fig/plot_SFx3.R 2_analysis/fst/50k/ \
  2_analysis/summaries/fst_outliers_998.tsv \
  2_analysis/summaries/fst_globals.txt

Rscript --vanilla R/fig/plot_SFx4.R 2_analysis/pi/50k/

Rscript --vanilla R/fig/plot_SFx5.R 2_analysis/pi/50k/ \
  2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz

Rscript --vanilla R/fig/plot_SFx6.R 2_analysis/dxy/50k/

Rscript --vanilla R/fig/plot_SFx7.R 2_analysis/fasttree/ no_outgroups.all.whg.SNP.tree