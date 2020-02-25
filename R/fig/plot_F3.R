#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F3.R \
#   2_analysis/msmc/output/ 2_analysis/cross_coalescence/output/ \
#   2_analysis/msmc/setup/msmc_grouping.txt 2_analysis/msmc/setup/msmc_cc_grouping.txt \
#   2_analysis/summaries/fst_globals.txt
# ===============================================================
# This script produces Figure 3 of the study "The genomic onset of a marine radiation"
# by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/msmc/output/', '2_analysis/cross_coalescence/output/',
# '2_analysis/msmc/setup/msmc_grouping.txt', '2_analysis/msmc/setup/msmc_cc_grouping.txt',
# '2_analysis/summaries/fst_globals.txt')
# ----------------------------------------
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(patchwork)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
msmc_path <- as.character(args[1])
cc_path <- as.character(args[2])
msmc_group_file <- as.character(args[3])
cc_group_file <- as.character(args[4])
fst_globals_file <- as.character(args[5])

# actual script =========================================================
msmc_groups <- read_tsv(msmc_group_file)
cc_groups <- read_tsv(cc_group_file)
fst_globals <- vroom::vroom(fst_globals_file,delim = '\t',
                            col_names = c('loc','run_prep','mean_fst','weighted_fst')) %>%
  separate(run_prep,into = c('pop1','pop2'),sep = '-') %>%
  mutate(run = str_c(pop1,loc,'-',pop2,loc),
         run = fct_reorder(run,weighted_fst))

msmc_files <- dir(msmc_path, pattern = '.final.txt.gz')
cc_files <- dir(cc_path, pattern = '.final.txt.gz')

msmc_data <- msmc_files %>%
  map_dfr(.f = get_msmc, msmc_path = msmc_path)

cc_data <- cc_files %>%
  map_dfr(get_cc, cc_groups = cc_groups, cc_path = cc_path) %>%
  mutate( run = factor(run, levels = levels(fst_globals$run)))

clr_alt <- clr
clr_alt['uni'] <- rgb(.8,.8,.8)
clr_ticks <- 'lightgray'

p_msmc <- msmc_data %>%
  # remove the two first and last time segments
  filter(!time_index %in% c(0:2,29:31)) %>%
  ggplot( aes(x=YBP, y=Ne, group = run_nr, colour = spec)) +
  # add guides for the logarithmic axes
  annotation_logticks(sides="tl", color = clr_ticks) +
  # add the msmc data as lines
  geom_line()+
  # set the color scheme
  scale_color_manual('Species',
                     values = clr_alt, label = sp_labs) +
  # format the x axis
  scale_x_log10(expand = c(0,0),
                breaks = c(10^3, 10^4, 10^5),
                position = 'top',
                labels = c("1-3 kya", "10-30 kya", "100-300 kya"),
                name = "Years Before Present") +
  # format the y axis
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                breaks = c(10^3,10^4,10^5,10^6)) +
  # format the color legend
  guides(colour = guide_legend(title.position = "top",
                               override.aes = list(alpha = 1, size=1),
                               nrow = 2, byrow=TRUE)) +
  # set the axis titles
  labs(x="Generations Before Present",
       y = expression(Effective~Population~Size~(italic(N[e])))) +
  # set plot range
  coord_cartesian(xlim = c(250, 5*10^5)) +
  # tune plot appreance
  theme_minimal()+
  theme(axis.ticks = element_line(colour = clr_ticks),
        legend.position = c(1,.03),
        legend.justification = c(1,0),
        legend.text.align = 0,
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        title = element_text(face = 'bold'),
        axis.title = element_text(face = 'plain'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(face = 'plain'))

p_cc <- cc_data %>%
  # remove the two first and last time segments
  filter( !time_index %in% c(0:2,29:31)) %>%
  arrange(run_nr) %>%
  # attach fst data
  left_join(fst_globals %>%
              select(run, weighted_fst)) %>%
  ggplot(aes(x = YBP, y = Cross_coal, group = run_nr, color = weighted_fst)) +
  # add guides for the logarithmic axis
  annotation_logticks(sides="b", color = clr_ticks) +
  # add the msmc data as lines
  geom_line(alpha = 0.2)+
  # set the color scheme
  scale_color_gradientn(name = expression(Global~weighted~italic(F[ST])),
                        colours = hypogen::hypo_clr_LGs[1:24])+
  # format the x axis
  scale_x_log10(expand = c(0,0),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # set the axis titles
  guides(color = guide_colorbar(barheight = unit(7, 'pt'),
                                barwidth = unit(150, 'pt'),
                                title.position = 'top'))+
  # set the axis titles
  labs(x = "Generations Before Present",
       y = 'Cross-coalescence Rate') +
  # set plot range
  coord_cartesian(xlim = c(250, 5*10^5)) +
  # tune plot appreance
  theme_minimal()+
  theme(axis.ticks = element_line(colour = clr_ticks),
        legend.position = c(1,.03),
        legend.direction = 'horizontal',
        legend.justification = c(1,0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        title = element_text(face = 'bold'),
        axis.title = element_text(face = 'plain'),
        legend.title = element_text(face = 'plain'))


p_done <- p_msmc / p_cc + plot_annotation(tag_levels = c('a'))


hypo_save(plot = p_done, filename = 'figures/F3.pdf',
          width = 10, height = 7,
          comment = 'Rscript srews this one up',
          device = cairo_pdf)
