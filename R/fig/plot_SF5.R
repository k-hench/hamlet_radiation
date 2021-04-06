#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF5.R 2_analysis/pi/50k/ \
#   2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
# ===============================================================
# This script produces Suppl. Figure 5 of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/pi/50k/',
#           '2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz')
# script_name <- "R/fig/plot_SF5.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(vroom)
library(hypoimg)
library(hypogen)
cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
pi_path <- as.character(args[1])
rho_path <- as.character(args[2])

# locate pi data files
files <- dir(pi_path, pattern = '^[a-z]{6}.50k')

# load pi data
data <- str_c(pi_path, files) %>%
  purrr::map(get_pi) %>%
  bind_rows()

# compute genome wide average pi for the subplot order
global_bar <- data %>%
  filter( BIN_START %% 50000 == 1) %>%
  select(N_VARIANTS, PI, spec) %>%
  group_by(spec) %>%
  summarise(genome_wide_pi = sum(N_VARIANTS*PI)/sum(N_VARIANTS)) %>%
  arrange(genome_wide_pi) %>%
  ungroup() %>%
  mutate(spec = fct_reorder(.f = spec, .x = genome_wide_pi),
         scaled_pi = genome_wide_pi/max(genome_wide_pi))

# load recombination data
rho_data <- vroom(rho_path, delim = '\t') %>%
  select(-BIN_END)

# merge pi and recombination data
combined_data <- data %>%
  # filter pi data to "non-overlapping" windows
  filter(BIN_START %% 50000 == 1 ) %>%
  # reorder populations by genome wide average pi
  mutate(spec = factor(spec, levels = levels(global_bar$spec))) %>%
  # merge with recombination data
  left_join(rho_data, by = c(CHROM = 'CHROM', BIN_START = 'BIN_START'))

# nest data to run linear regression on all runs in one go
model_data <- combined_data %>%
  group_by(spec) %>%
  nest() %>%
  left_join(., global_bar) %>%
  mutate(mod =  map(data, ~ lm(.$PI ~ .$RHO))) %>%
  bind_cols(., summarise_model(.))

# create table with fish annotations
grob_tibble2 <- global_bar$spec %>%
  purrr::map(fish_plot2) %>%
  bind_rows()

# compose final figure
p <- combined_data %>%
  ggplot()+
  # add fish annotations
  geom_hypo_grob2(data = grob_tibble2,
                  aes(grob = grob, rel_x = .25,rel_y = .75),
                  angle = 0, height = .5,width = .5)+
  # add hex-bin desity layer
  geom_hex(bins = 30,color = rgb(0,0,0,.3),
           aes(fill=log10(..count..), x = RHO, y = PI))+
 # general plot structure (separated by run)
  facet_wrap(spec ~., ncol = 3)+
  # set axis layout and color scheme
  scale_x_continuous(name = expression(rho))+
  scale_y_continuous(name = expression(pi))+
  scico::scale_fill_scico(palette = 'berlin') +
  # customize legend
  guides(fill = guide_colorbar(direction = 'horizontal',
                               title.position = 'top',
                               barheight = unit(7,'pt'),
                               barwidth = unit(130,'pt')))+
  # general plot layout
  theme_minimal()+
  theme(legend.position = c(.84,.01),
        strip.text = element_blank())

# export final figure
hypo_save(filename = 'figures/SF5.pdf',
          plot = p,
          width = 8,
          height = 10,
          comment = plot_comment)
