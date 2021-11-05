#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF21.R \
#     2_analysis/pi/50k/
# ===============================================================
# This script produces Suppl. Figure 21 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/pi/50k/')
# script_name <- "R/fig/plot_SF21.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(hypoimg)
library(hypogen)
library(ggtext)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
pi_path <- as.character(args[1])

# locate pi data files
files <- dir(pi_path, pattern = '^pi.[a-z]{6}.50k')
# load pi data
data <- str_c(pi_path, files) %>%
  purrr::map(get_pi) %>%
  bind_rows()

# create table for the indication of genome wide average pi in the plot background
# (rescale covered pi range to the extent of the genome)
global_bar <- data %>%
  # filter to non-overlaping windows only
  filter( BIN_START %% 50000 == 1) %>%
  select(N_SITES, PI, spec) %>%
  group_by(spec) %>%
  summarise(genome_wide_pi = sum(N_SITES*PI)/sum(N_SITES)) %>%
  arrange(genome_wide_pi) %>%
  ungroup() %>%
  mutate(spec = fct_reorder(.f = spec, .x = genome_wide_pi),
         scaled_pi = genome_wide_pi/max(genome_wide_pi))

global_bar %>% write_tsv("2_analysis/summaries/pi_globals.tsv")

# prepare plot annotaton images
grob_tibble <- global_bar$spec %>%
  purrr::map(fish_plot) %>%
  bind_rows()

# prepare plotting elements --------
# pre-define secondary x-axis breaks
sc_ax <- scales::cbreaks(c(0,max(global_bar$genome_wide_pi)),
                         scales::pretty_breaks(4))

# pre-define secondary x-axis labels
labels <- str_c(c("", sc_ax$breaks[2:6]*1000),
                c("0", rep("\u00B710^-3",5)))

# sort pair-wise population comparisons by average genome wide pi
data <- data %>%
  mutate(spec = factor(spec, levels = levels(global_bar$spec)))

# compose final figure
p_done <- ggplot()+
  # general plot structure separated by run
  facet_wrap( .~spec, as.table = TRUE, ncol = 1, dir = 'v')+
  # add genome wide average pi in the background
  geom_rect(data = global_bar %>% mutate(xmax = scaled_pi * hypo_karyotype$GEND[24]),
            aes(xmin = 0, xmax = xmax, ymin = -Inf, ymax = Inf),
            color = rgb(1,1,1,0),
            fill = clr_below)+
  # add LG borders
  geom_vline(data = hypogen::hypo_karyotype,aes(xintercept = GEND),color = hypo_clr_lg)+
  # add pi data points
  geom_point(data = data, aes(x = gpos, y = PI),
             size=.2, color = plot_clr) +
  # add fish images
  geom_hypo_grob2(data = grob_tibble,
                  aes(grob = grob, rel_x = .975, rel_y = .5),
                  angle = 0, height = .8, width = .12)+
  # set axis layout
  scale_x_hypo_LG(sec.axis =  sec_axis(~ ./hypo_karyotype$GEND[24],
                                       breaks = (sc_ax$breaks/max(global_bar$genome_wide_pi)),
                                       labels = labels,
                                       name = "Genomic position/ Genome wide *\u03C0*"))+
  scale_y_continuous(name = "*\u03C0*", breaks = c(0,.01,.02))+
  # set plot extent
  coord_cartesian(xlim = c(0, hypo_karyotype$GEND[24]*1.06))+
  # general plot layout
  theme_hypo()+
  theme(strip.text = element_blank(),
        legend.position = 'none',
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        axis.text.x.bottom = element_markdown(colour = 'darkgray'))

# export final figure
hypo_save(filename = 'figures/SF21.png',
          plot = p_done,
          width = 8,
          height = 8,
          dpi = 600,
          type = "cairo",
          comment = plot_comment)

system("convert figures/SF21.png figures/SF21.pdf")
system("rm figures/SF21.png")
create_metadata <- str_c("exiftool -overwrite_original -Description=\"", plot_comment, "\" figures/SF21.pdf")
system(create_metadata)

# export S.Tab. 4
global_bar %>% 
  dplyr::select(- scaled_pi) %>% 
  mutate(loc = loc_names[str_sub(spec, -3,-1)],
         spec = str_c("H. ", sp_names[str_sub(spec, 1,3)]),
         genome_wide_pi = sprintf("%.5f",genome_wide_pi)) %>% 
  pivot_wider(names_from = spec, values_from = genome_wide_pi, values_fill = "-") %>% 
  knitr::kable(format = "latex")
