#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F3.R \
#    2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt
# ===============================================================
# This script produces Figure 3 of the study "Ancestral variation, hybridization and modularity
# fuel a marine radiation" by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/fst/50k/', '2_analysis/summaries/fst_globals.txt')
# script_name <- "plot_F3.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(ggforce)
library(hypoimg)
library(vroom)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)
# config -----------------------
data_dir <- as.character(args[1])
globals_file <- as.character(args[2])
# script -----------------------

# locate data files
files <- dir(path = data_dir, pattern = '.50k.windowed.weir.fst.gz')

# load genome wide average fst data
globals <- vroom::vroom(globals_file, delim = '\t',
                        col_names = c('loc','run','mean','weighted')) %>%
  mutate(run = str_c(loc,'-',run) %>%
           reformat_run_name()
  )

# prepare data import settings within a data table (tibble)
import_table <- list(file = str_c(data_dir,files),
                     fst_threshold = c(.5,.4,.3,.2,.1,.05,.02,.01)) %>%
  cross_df() %>%
  mutate( run =  file %>%
            str_remove('^.*/') %>%
            str_sub(.,1,11) %>%
            reformat_run_name())

# load data and compute statistics based on fixed fst treshold
data <- purrr::pmap_dfr(import_table,get_fst_fixed) %>%
  left_join(globals) %>%
  mutate(run = fct_reorder(run, weighted))


# pre-format labels
data2 <- data %>%
  select(threshold_value,weighted,n,avg_length,overal_length) %>%
  mutate(avg_length = avg_length/1000,
         overal_length = overal_length/(10^6)) %>%
  rename(`Number~of~Regions` = 'n',
         `Average~Region~Length~(kb)` = 'avg_length',
         `Cummulative~Region~Length~(Mb)` = 'overal_length') %>%
  pivot_longer(names_to = 'variable',values_to = 'Value',3:5) %>%
  mutate(threshold_value = str_c('italic(F[ST])~threshold:~',
                                 threshold_value),
         variable = factor(variable, levels = c('Number~of~Regions',
                                                'Average~Region~Length~(kb)',
                                                'Cummulative~Region~Length~(Mb)')))

# set font size
fnt_sz <- 10

# compiule plot
p <-  data2 %>%
  # select thresholds of interest
  filter(!(threshold_value %in% (c(0.02,.1,0.3, .4) %>%
                                   str_c("italic(F[ST])~threshold:~",.)))) %>%
  ggplot(aes(x = weighted, y = Value, fill = weighted))+
  # add red line for genome extent in lowest row
  geom_hline(data = tibble(variable = factor(c('Cummulative~Region~Length~(Mb)',
                                               'Average~Region~Length~(kb)',
                                               'Number~of~Regions'),
                                             levels = c('Number~of~Regions',
                                                        'Average~Region~Length~(kb)',
                                                        'Cummulative~Region~Length~(Mb)')),
                           y = c(559649677/(10^6),NA,NA)),
             aes(yintercept = y),
             color=rgb(1,0,0,.25))+
  # add data points
  geom_point(size = 1.75,
             color = plot_clr, shape = 21)+
  # define plot stucture
  facet_grid(variable~threshold_value,
             scale='free',
             switch = 'y',
             labeller = label_parsed)+
  # configure scales
  scale_fill_gradientn(name = expression(weighted~italic(F[ST])),
                       colours = hypogen::hypo_clr_LGs[1:24] %>% clr_lighten(factor = .3))+
  scale_x_continuous(name = expression(Whole-genome~differentiation~(weighted~italic(F[ST]))),
                     breaks = c(0,.05,.1),
                     limits = c(-.00025,.10025),
                     labels = c("0", "0.05", "0.1"))+
  # configure legend
  guides(fill = guide_colorbar(barwidth = unit(150, "pt"),
                               label.position = "top",
                               barheight = unit(5,"pt")))+
  # tweak plot apperance
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = .5, angle = 0),
        axis.title.x = element_text(vjust = -2),
        legend.position = "bottom",
        strip.text = element_text(size = fnt_sz),
        legend.direction = "horizontal",
        strip.placement = 'outside',
        axis.title = element_text(size = fnt_sz),
        legend.title = element_text(size = fnt_sz),
        strip.background.y = element_blank())

# export figure 3
scl <- .675
hypo_save(filename = 'figures/F3.pdf',
          plot = p,
          width = 16*scl,
          height = 12*scl,
          device = cairo_pdf,
          comment = plot_comment)
