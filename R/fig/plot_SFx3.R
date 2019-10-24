#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFx3.R figures/data/fst/ \
#   figures/data/summaries/all_multi_fst_outliers_998.tsv \
#   figures/data/summaries/fst_globals.txt
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c('figures/data/fst/', 'figures/data/summaries/all_multi_fst_outliers_998.tsv', 'figures/data/summaries/fst_globals.txt')
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(vroom)

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

data <- purrr::pmap(tibble(file = str_c(data_path,files),run = run_files),
                    hypo_import_windows) %>%
  bind_rows() %>%
  set_names(., nm = c('CHROM', 'BIN_START', 'BIN_END', 'N_VARIANTS', 'WEIGHTED_FST', 'MEAN_FST', 'GSTART', 'POS', 'GPOS', 'run'))

# global bars ------------------------
global_bar <- globals %>%
  select(weighted,run) %>%
  mutate(run = as.character(run)) %>%
  setNames(.,nm = c('fst','run')) %>%
  pmap(.,fst_bar_row_run) %>%
  bind_rows()

# annotaton ---------------------------
#### -------------- vvvvvvv -----------
runs <- globals %>%
  group_by(run) %>%
  count() %>%
  ungroup() %>%
  select(-n) %>%
  mutate(loc = str_sub(run,4,6),
         right_short = str_sub(run,1,3),
         left_short = str_sub(run,8,10)) %>%
  mutate(left = left_short,
         right = right_short)

grob_tibble <- runs %>%
  select(1:2, 5:6) %>%
  pmap(.,plot_pair_run) %>%
  bind_rows()

# arrange factor order ----------------
data$run <- refactor_run(data, globals)
grob_tibble$run <- refactor_run(grob_tibble, globals)
global_bar$run <- refactor_run(global_bar,globals)

# plotting ---------------
sc_ax <- scales::cbreaks(c(0,max(globals$weighted)),
                         scales::pretty_breaks(4))

p <- ggplot()+
  facet_wrap( . ~ run, as.table = TRUE, ncol = 2,dir = 'v') +
  geom_rect(data = global_bar %>%
              mutate(xmax = xmax * hypo_karyotype$GEND[24]),
            aes(xmin = 0, xmax = xmax,
                ymin = -Inf, ymax = Inf),
            color = rgb(1,1,1,0),
            fill = clr_below) +
  geom_vline(data = hypogen::hypo_karyotype,
             aes(xintercept = GEND),
             color = hypo_clr_lg) +
  geom_hypo_grob(data = grob_tibble,
                 aes(grob = grob, x = .9,y = .7),
                 angle = 0, height = .5, width = .16) +
  geom_point(data = data, aes(x = GPOS, y = WEIGHTED_FST),
             size=.2,color = plot_clr) +
  scale_x_hypo_LG(sec.axis =  sec_axis(~ ./hypo_karyotype$GEND[24],
                                       breaks = (sc_ax$breaks/max(globals$weighted)),
                                       labels = sprintf("%.2f", sc_ax$breaks),
                                       name = expression(Genomic~position/~Genome~wide~weighted~italic(F[ST])))) +
  scale_y_continuous(name = expression(italic('F'[ST])),
                     limits = c(-.1,1)) +
  theme_hypo() +
  theme(strip.text = element_blank(),
        legend.position = 'none',
        axis.title.x = element_text(),
        axis.text.x.bottom = element_text(colour = 'darkgray'))

hypo_save(filename = 'figures/SX3.png',
          plot = p,
          width = 16,
          height = 9,
          comment = plot_comment)
