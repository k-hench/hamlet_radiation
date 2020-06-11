#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF8.R 2_analysis/pi/50k/
# ===============================================================
# This script produces Suppl. Figure 8 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/pi/50k/')
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
pi_path <- as.character(args[1])

# load data -------------------
files <- dir(pi_path, pattern = '^[a-z]{6}.50k')
data <- str_c(pi_path, files) %>%
  purrr::map(get_pi) %>%
  bind_rows()

# global bars ------------------------
global_bar <- data %>%
  filter( BIN_START %% 50000 == 1) %>%
  select(N_VARIANTS, PI, spec) %>%
  group_by(spec) %>%
  summarise(genome_wide_pi = sum(N_VARIANTS*PI)/sum(N_VARIANTS)) %>%
  arrange(genome_wide_pi) %>%
  ungroup() %>%
  mutate(spec = fct_reorder(.f = spec, .x = genome_wide_pi),
         scaled_pi = genome_wide_pi/max(genome_wide_pi))

global_bar %>%
  summarise(avg_pi = mean(genome_wide_pi),
            sd_pi = sd(genome_wide_pi))
# annotaton ---------------------------
grob_tibble <- global_bar$spec %>%
  purrr::map(fish_plot) %>%
  bind_rows()

# plotting ---------------
sc_ax <- scales::cbreaks(c(0,max(global_bar$genome_wide_pi)), scales::pretty_breaks(4))

labels <-  c("0",
             substitute(a%.%10^-3, list(a = sprintf("%.0f", sc_ax$breaks[2]*1000))),
             substitute(a%.%10^-3, list(a = sprintf("%.0f", sc_ax$breaks[3]*1000))),
             substitute(a%.%10^-3, list(a = sprintf("%.0f", sc_ax$breaks[4]*1000))),
             substitute(a%.%10^-3, list(a = sprintf("%.0f", sc_ax$breaks[5]*1000))),
             substitute(a%.%10^-3, list(a = sprintf("%.0f", sc_ax$breaks[6]*1000))))

data <- data %>%
  mutate(spec = factor(spec, levels = levels(global_bar$spec)))

p <- ggplot()+
  facet_wrap( .~spec, as.table = TRUE, ncol = 1,dir = 'v')+
  geom_rect(data = global_bar %>% mutate(xmax = scaled_pi * hypo_karyotype$GEND[24]),
            aes(xmin = 0, xmax = xmax, ymin = -Inf, ymax = Inf),
            color = rgb(1,1,1,0),
            fill = clr_below)+
  geom_vline(data = hypogen::hypo_karyotype,aes(xintercept = GEND),color = hypo_clr_lg)+
  geom_point(data = data, aes(x = gpos, y = PI),
             size=.2,color = plot_clr) +
  geom_hypo_grob2(data = grob_tibble,
                  aes(grob = grob, rel_x = .975,rel_y = .5),
                  angle = 0, height = .8,width = .12)+
  scale_x_hypo_LG(sec.axis =  sec_axis(~ ./hypo_karyotype$GEND[24],
                                       breaks = (sc_ax$breaks/max(global_bar$genome_wide_pi)),
                                       labels = labels,
                                       name = expression(Genomic~position/~Genome~wide~italic(pi))))+
  scale_y_continuous(name = expression(italic(pi)), breaks = c(0,.006,.012))+
  coord_cartesian(xlim = c(0,hypo_karyotype$GEND[24]*1.06))+
  theme_hypo()+
  theme(strip.text = element_blank(),
        legend.position = 'none',
        axis.title.x = element_text(),
        axis.text.x.bottom = element_text(colour = 'darkgray'))

hypo_save(filename = 'figures/SF8.png',
          plot = p,
          width = 8,
          height = 8,
          comment = plot_comment)
