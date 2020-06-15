#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF2.R 2_analysis/fst/50k/ \
#   2_analysis/summaries/fst_outliers_998.tsv \
#   2_analysis/summaries/fst_globals.txt
# ===============================================================
# This script produces Suppl. Figure 2 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/fst/50k/',
#           '2_analysis/summaries/fst_outliers_998.tsv',
#           '2_analysis/summaries/fst_globals.txt')
# script_name <- "R/fig/plot_SF2.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(vroom)
library(ggtext)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
data_path <- as.character(args[1])
outlier_file <- as.character(args[2])
globals_file <- as.character(args[3])
# load data -------------------
files <- dir(data_path,pattern = '.50k.windowed.weir.fst.gz')
run_files <- files %>%
  str_sub(.,1,11) %>%
  str_replace(.,pattern = '([a-z]{3})-([a-z]{3})-([a-z]{3})', '\\2\\1-\\3\\1')

globals <- vroom::vroom(globals_file, delim = '\t',
                        col_names = c('loc','run','mean','weighted')) %>%
  separate(run, into = c('pop1','pop2')) %>%
  mutate(run = str_c(pop1,loc,'-',pop2,loc),
         run = fct_reorder(run,weighted))

data <- purrr::pmap(tibble(file = str_c(data_path,files),
                           run = run_files),
                    hypo_import_windows) %>%
  bind_rows() %>%
  set_names(., nm = c('CHROM', 'BIN_START', 'BIN_END', 'N_VARIANTS',
                      'WEIGHTED_FST', 'MEAN_FST', 'GSTART', 'POS', 'GPOS', 'run')) %>%
  mutate(pop1 = str_sub(run,1,3),
         pop2 = str_sub(run,8,10),
         loc = str_sub(run,4,6),
         run_label = str_c("*H. ", sp_names[pop1],"* - *H. ", sp_names[pop2],"*<br>(",loc_names[loc],")" ))

# global bars ------------------------
global_bar <- globals %>%
  select(weighted,run) %>%
  mutate(run = as.character(run)) %>%
  setNames(.,nm = c('fst','run')) %>%
  pmap(.,fst_bar_row_run) %>%
  bind_rows() %>%
  mutate(pop1 = str_sub(run,1,3),
         pop2 = str_sub(run,8,10),
         loc = str_sub(run,4,6),
         run_label = str_c("*H. ", sp_names[pop1],"* - *H. ", sp_names[pop2],"*<br>(",loc_names[loc],")" ),
         run_label = fct_reorder(run_label,xmax_org))

# plotting ---------------
sc_ax <- scales::cbreaks(c(0,max(globals$weighted)),
                         scales::pretty_breaks(4))

p <- ggplot()+
  facet_grid( run_label ~ ., as.table = TRUE) +
  geom_rect(data = global_bar %>%
              mutate(xmax = xmax * hypo_karyotype$GEND[24]),
            aes(xmin = 0, xmax = xmax,
                ymin = -Inf, ymax = Inf),
            color = rgb(1,1,1,0),
            fill = clr_below) +
  geom_vline(data = hypogen::hypo_karyotype,
             aes(xintercept = GEND),
             color = hypo_clr_lg) +
  geom_point(data = data  %>% mutate(run_label = factor(run_label, levels = levels(global_bar$run_label))),
              aes(x = GPOS, y = WEIGHTED_FST),
              size=.2,color = plot_clr) +
  scale_x_hypo_LG(sec.axis =  sec_axis(~ ./hypo_karyotype$GEND[24],
                                       breaks = (sc_ax$breaks/max(globals$weighted)),
                                       labels = sprintf("%.2f", sc_ax$breaks),
                                       name = expression(Genomic~position/~Genome~wide~weighted~italic(F[ST])))) +
  scale_y_continuous(name = expression(italic('F'[ST])),
                     limits = c(-.1,1),
                     breaks = c(0,.5,1)) +
  theme_hypo() +
  theme(strip.text.y = element_markdown(angle = 0),
        strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_text(),
        axis.text.x.bottom = element_text(colour = 'darkgray'))

hypo_save(filename = 'figures/SF2.png',
          plot = p,
          width = 8,
          height = 12,
          comment = plot_comment)
