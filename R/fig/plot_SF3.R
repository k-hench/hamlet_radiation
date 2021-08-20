#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF3.R \
#     2_analysis/pca/
# ===============================================================
# This script produces Suppl. Figure 3 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/pca/")
# script_name <- "R/fig/plot_SF3.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(patchwork)
library(hypoimg)
library(hypogen)
cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)

# config -----------------------
pca_dir <- as.character(args[1])
clr_alt <- clr
clr_alt[["uni"]] <- "lightgray"

fish_tib <- tibble(short = names(clr)[!names(clr) %in% c("flo", "tab", "tor")],
                   x = c(0.5,  3.5,  7,  9.7, 12.25, 15.25, 18, 21.5))

key_sz <- .75
sp_fam <- rep(c("H", "S", "H"), c(8, 2, 1)) %>% set_names(nm = names(sp_names))
p_leg <- fish_tib %>% 
  ggplot() +
  coord_equal(xlim = c(-.05, 24), expand = 0) +
 geom_tile(aes(x = x, y = 0,
                fill = short, 
                color = after_scale(prismatic::clr_darken(fill, .25))),
            width = key_sz, height = key_sz, size = .3) +
  geom_text(aes(x = x + .6, y = 0,
                label = str_c(sp_fam[short], ". ", sp_names[short])), 
            hjust = 0, fontface = "italic", size = plot_text_size / ggplot2:::.pt) +
  pmap(fish_tib, plot_fish_lwd, width = 1, height = 1, y = 0) +
  scale_fill_manual(values = clr, guide = FALSE) +
  theme_void()

p_done <- cowplot::plot_grid((tibble(loc = c("bel.", "hon.", "pan."), 
                              mode = rep(c("subset_non_diverged"), 3),
                              pc1 = 1,
                              pc2 = 2) %>% 
                        pmap(pca_plot_no_fish) %>% 
                        wrap_plots(ncol = 3) +
                        plot_annotation(tag_levels = "a") & 
                        theme(plot.background = element_blank())),
                     p_leg,
  ncol = 1 ,
  rel_heights = c(1,.1))

hypo_save(p_done, filename = 'figures/SF3.pdf',
          width = f_width,
          height = f_width * .38,
          device = cairo_pdf,
          bg = "transparent",
          comment = plot_comment)
