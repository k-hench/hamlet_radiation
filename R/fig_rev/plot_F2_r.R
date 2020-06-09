#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_F2.R 2_analysis/fst/50k/ 2_analysis/summaries/fst_globals.txt
# ===============================================================
# This script produces Figure 1 of the study "The genomic onset of a marine radiation"
# by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/fst/50k/', '2_analysis/summaries/fst_globals.txt')
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

# load data -------------------
files <- dir(path = data_dir, pattern = '.50k.windowed.weir.fst.gz')

globals <- vroom::vroom(globals_file, delim = '\t',
                        col_names = c('loc','run','mean','weighted')) %>%
  mutate(run = str_c(loc,'-',run) %>%
           reformat_run_name()#,
         #run = fct_reorder(run,weighted)
  )

# fixed fst -----------
import_table <- list(file = str_c(data_dir,files),
                     fst_threshold = c(.6,.4,.3,.2,.1,.05,.02,.01)) %>%
  #                   fst_threshold = c(.75,.7,.65,.6,.55,.5,.45,.4,.35,.3,.25,.2)) %>%
  cross_df() %>%
  mutate( run =  file %>%
            str_remove('^.*/') %>%
            str_sub(.,1,11) %>%
            reformat_run_name())

data <- purrr::pmap_dfr(import_table,get_fst_fixed) %>%
  left_join(globals) %>%
  mutate(run = fct_reorder(run, weighted))

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

fnt_sz <- 10
p <-  data2 %>%
 filter(!(threshold_value %in% (c(0.02,.1,0.3, .6) %>%
        str_c("italic(F[ST])~threshold:~",.)))) %>%
  ggplot(aes(x = weighted, y = Value, fill = weighted))+
  geom_hline(data = tibble(variable = factor(c('Cummulative~Region~Length~(Mb)',
                                               'Average~Region~Length~(kb)',
                                               'Number~of~Regions'), levels = c('Number~of~Regions',
                                                                                'Average~Region~Length~(kb)',
                                                                                'Cummulative~Region~Length~(Mb)')),
                           y = c(559649677/(10^6),NA,NA)),
             aes(yintercept = y),
             color=rgb(1,0,0,.25))+
  geom_point(size = 1.75,#.8,
             color = plot_clr, shape = 21)+
  facet_grid(variable~threshold_value,
             scale='free',#nrow = 3,
             switch = 'y',
             labeller = label_parsed)+
  scale_fill_gradientn(name = expression(weighted~italic(F[ST])),
                       colours = hypogen::hypo_clr_LGs[1:24] %>% clr_lighten(factor = .3))+
  # scale_fill_viridis_c(name = expression(weighted~italic(F[ST])),
  #                      option = "B")+
  scale_x_continuous(name = expression(Whole-genome~differentiation~(weighted~italic(F[ST]))),
                     breaks = c(0,.05,.1),
                     limits = c(-.00025,.10025),
                     labels = c("0", "0.05", "0.1"))+
#  theme_minimal()+
  guides(fill = guide_colorbar(barwidth = unit(150, "pt"),
                               label.position = "top",
                                barheight = unit(5,"pt")))+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = .5, angle = 0),
        legend.position = c(1,1),
        legend.background = element_rect(color = plot_clr,fill = rgb(1,1,1,.3), size = .2),
        legend.justification = c(1,1),
        strip.text = element_text(size = fnt_sz),
        legend.direction = "horizontal",
        strip.placement = 'outside',
        axis.title = element_text(size = fnt_sz),
        legend.title = element_text(size = fnt_sz),
     #   panel.spacing.x = unit(10,"pt"),
        strip.background.y = element_blank())

scl <- .675
hypo_save(filename = 'figures/F2.pdf',
          plot = p,
          width = 16*scl,
          height = 11*scl,
          device = cairo_pdf,
          comment = plot_comment)

