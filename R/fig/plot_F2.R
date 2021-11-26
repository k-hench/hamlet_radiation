#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F2.R \
#     2_analysis/fst/50k/multi_fst.50k.tsv.gz \
#     2_analysis/GxP/50000/ \
#     2_analysis/summaries/fst_outliers_998.tsv \
#     https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R \
#     2_analysis/twisst/weights/ \
#     ressources/plugin/trees/ \
#     2_analysis/summaries/fst_globals.txt
# ===============================================================
# This script produces Figure 2 of the study "Rapid radiation in a highly
# diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/fst/50k/multi_fst.50k.tsv.gz',
# '2_analysis/GxP/50000/', '2_analysis/summaries/fst_outliers_998.tsv',
# 'https://raw.githubusercontent.com/simonhmartin/twisst/master/plot_twisst.R',
# '2_analysis/twisst/weights/', 'ressources/plugin/trees/',
# '2_analysis/summaries/fst_globals.txt')
# script_name <- "R/fig/plot_F2.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
fst_file <- as.character(args[1])
gxp_dir <- as.character(args[2])
outlier_table <- as.character(args[3])
twisst_script <- as.character(args[4])
w_path <- as.character(args[5])
d_path <- as.character(args[6])
global_fst_file <- as.character(args[7])
source(twisst_script)

# start script -------------------
# import fst data
fst_data <- vroom::vroom(fst_file, delim = '\t') %>%
  select(CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST) %>%
  setNames(., nm = c('CHROM', 'BIN_START', 'BIN_END', 'n_snps', 'fst') ) %>%
  add_gpos() %>%
  select(GPOS, fst) %>%
  setNames(., nm = c('GPOS','value')) %>%
  mutate(window = str_c('bold(',project_case('a'),'):joint~italic(F[ST])'))

# set G x P traits to be imported
traits <- c("Bars.lm.50k.5k.txt.gz", "Peduncle.lm.50k.5k.txt.gz", "Snout.lm.50k.5k.txt.gz")

# set trait figure panels
trait_panels <- c(Bars = str_c('bold(',project_case('d'),')'),
                  Peduncle = str_c('bold(',project_case('e'),')'),
                  Snout = str_c('bold(',project_case('f'),')'))

# import G x P data
gxp_data <- str_c(gxp_dir,traits) %>%
  purrr::map(get_gxp) %>%
  join_list() %>%
  gather(key = 'window', value = 'value',2:4)

# import genome wide Fst data summary  --------
globals <- vroom::vroom(global_fst_file, delim = '\t',
                        col_names = c('loc','run','mean','weighted')) %>%
  mutate(run = str_c(str_sub(run,1,3),loc,'-',str_sub(run,5,7),loc),
         run = fct_reorder(run,weighted))

# dxy and pi are only shown for one exemplary population (/pair)
# select dxy pair run (15 is one of the two central runs of the 28 pairs)
# here, the 15th lowest fst value is identified as "selector"
selectors_dxy <- globals %>%
  arrange(weighted) %>%
  .$weighted %>%
  .[15]


# import topology weighting data
twisst_data <- tibble(loc = c('bel','hon'),
                      panel = c('b','c') %>% project_case() %>% str_c('bold(',.,')')) %>%
  purrr::pmap(match_twisst_files) %>%
  bind_rows() %>%
  select(GPOS, topo3,topo_rel,window,weight)

# the "null-weighting" is computed for both locations
twisst_null <- tibble(window = c(str_c('bold(',project_case('b'),'):~italic(w)[bel]'),
                                 str_c('bold(',project_case('c'),'):~italic(w)[hon]')),
                      weight = c(1/15, 1/105))

# combine data types --------
data <- bind_rows(fst_data, gxp_data)

# import fst outliers
outliers <-  vroom::vroom(outlier_table, delim = '\t')

# the focal outlier IDs are set
outlier_pick <- c('LG04_1', 'LG12_3', 'LG12_4')

# the table for the outlier labels is created
outlier_label <- outliers %>%
  filter(gid %in% outlier_pick) %>%
  mutate(label = letters[row_number()] %>% project_inv_case(),
         x_shift_label = c(-1,-1.2,1)*10^7,
         gpos_label = gpos + x_shift_label,
         gpos_label2 = gpos_label - sign(x_shift_label) *.5*10^7,
         window = str_c('bold(',project_case('a'),'):joint~italic(F[ST])'))

# the y height of the outlier labels and the corresponding tags is set
outlier_y <- .45
outlier_yend <- .475

# the icons for the traits of the GxP are loaded
trait_tibble <- tibble(window = c("bold(d):italic(p)[Bars]",
                                  "bold(e):italic(p)[Peduncle]",
                                  "bold(f):italic(p)[Snout]"),
                       grob = hypo_trait_img$grob_circle[hypo_trait_img$trait %in% c('Bars', 'Peduncle', 'Snout')])

# finally, the figure is being put together
p_done <- ggplot()+
  # add gray/white LGs background
  geom_hypo_LG()+
  # the red highlights for the outlier regions are added
  geom_vline(data = outliers, aes(xintercept = gpos), color = outlr_clr)+
  # the tags of the outlier labels are added
  geom_segment(data = outlier_label,
               aes(x = gpos,
                   xend = gpos_label2, y = outlier_y, yend = outlier_yend),
               color = alpha(outlr_clr,1),
               size = .2)+
  # the outlier labels are added
  geom_text(data = outlier_label, aes(x = gpos_label, y = outlier_yend, label = label),
            color = alpha(outlr_clr,1), 
            fontface = 'bold',
            size = plot_text_size / ggplot2:::.pt)+
  # the fst, delta dxy and gxp data is plotted
  geom_point(data = data, aes(x = GPOS, y = value),size = plot_size, color = plot_clr) +
  # the topology weighting data is plotted
  geom_line(data = twisst_data, aes(x = GPOS, y = weight, color = topo_rel), size = .4) +
  # the null weighting is added
  geom_hline(data = twisst_null, aes(yintercept = weight), color = rgb(1, 1, 1, .5), size = .4) +
  # the trait icons are added
  geom_hypo_grob(data = trait_tibble,
                 aes(grob = grob, angle = 0, height = .65),
                 inherit.aes = FALSE, x = .95, y = 0.65)+
  # setting the scales
  scale_fill_hypo_LG_bg() +
  scale_x_hypo_LG()+
  scale_color_gradient( low = "#f0a830ff", high = "#084082ff", guide = FALSE)+
  # organizing the plot across panels
  facet_grid(window~.,scales = 'free',switch = 'y', labeller = label_parsed)+
  # tweak plot appreance
  theme_hypo()+
  theme(text = element_text(size = plot_text_size),
        legend.position = 'bottom',
        axis.title = element_blank(),
        strip.text = element_text(size = plot_text_size),
        strip.background = element_blank(),
        strip.placement = 'outside')

hypo_save(p_done,
          filename = 'figures/F2.png',
          width = f_width,
          height = f_width * .5,
          dpi = 600,
          type = "cairo",
          comment = plot_comment)

system("convert figures/F2.png figures/F2.pdf")
system("rm figures/F2.png")
create_metadata <- str_c("exiftool -overwrite_original -Description=\"", plot_comment, "\" figures/F2.pdf")
system(create_metadata)
