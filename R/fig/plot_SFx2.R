#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFx2.R \
#   metadata/phenotypes.sc ressources/img/
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c('metadata/phenotypes.sc', 'ressources/img/')
args = commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypoimg)
library(ggalluvial)

cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

cli::rule( left = str_c(crayon::bold('Script: '),crayon::red(script_name)))
args = args[7:length(args)]
cat(' ')
cat(str_c(crayon::green(cli::symbol$star),' ', 1:length(args),': ',crayon::green(args),'\n'))
cli::rule(right = getwd())

# config -----------------------
pheno_file <- as.character(args[1])
img_dir <- as.character(args[2])
wdh <- .15

bars <- hypo_read_svg(str_c(img_dir,'bars_r.c.svg'))
snout <- hypo_read_svg(str_c(img_dir,'snout_r.c.svg'))
peduncle <- hypo_read_svg(str_c(img_dir,'peduncle_r.c.svg'))

pheno_data <- read_sc(pheno_file) %>%
  select(spec, Bars, Peduncle, Snout, Lines,Tail_transparent, Blue ) %>%
  filter(!is.na(Bars)) %>%
  mutate(Bars = c('no', 'yes')[Bars+1],
         Peduncle = c('no', 'yes')[Peduncle+1],
         Snout = c('no', 'yes')[Snout+1],
         Lines = c('no', 'yes')[Lines+1],
         Tail_transparent = c('no', 'yes')[Tail_transparent+1],
         Blue = c('no', 'yes')[Blue+1]) %>%
  group_by(spec, Bars, Peduncle, Snout) %>% #, Lines,Tail_transparent, Blue) %>%
  count()

p_done <- pheno_data %>%
  ungroup() %>%
  ggplot(aes(axis1 = Bars,
             axis2 = Peduncle,
             axis3 = Snout,
             y = n)) +
  geom_alluvium(aes(color = spec,fill = spec,
                    group = spec),width = wdh, alpha = .65) +
  geom_stratum(width = wdh) +
  geom_text(stat = "stratum", label.strata = TRUE) +
  annotation_custom(grob = bars, xmin = .55,xmax = 1.5,ymin = 110,ymax = 150)+
  annotation_custom(grob = peduncle, xmin = 1.55,xmax = 2.5,ymin = 110,ymax = 150)+
  annotation_custom(grob = snout, xmin = 2.55,xmax = 3.5,ymin = 110,ymax = 150)+
  scale_x_discrete(limits = c("Bars", "Peduncle", "Snout"),#, 'Lines','Tail_transparent', 'Blue'),
                   expand = c(.095, .033),
                   position = 'top') +
  xlab("Trait") +
  scale_fill_manual('Species',
                    values = clr,
                    labels = sp_labs[names(clr)])+
  scale_color_manual('Species',
                     values = clr2,
                     labels = sp_labs[names(clr2)])+
  scale_y_continuous(expand = c(0, 0))+
  theme_minimal()+
  theme(legend.position = 'bottom',
        legend.text.align = 0,
        panel.grid = element_blank(),
        plot.background = element_blank())


hypo_save(p_done, filename = 'figures/SX2.pdf',
          width = 10, height = 5.5,
          comment = plot_comment)