#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFx3.R \
#     2_analysis/summaries/fst_outliers_998.tsv
# ===============================================================
# This script produces Suppl. Figure x1 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/summaries/fst_outliers_998.tsv")
# script_name <- "R/fig/plot_SFx1.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(tidygraph)
library(ggraph)
library(prismatic)
library(patchwork)
library(IRanges)
library(plyranges)
library(hypogen)
library(hypoimg)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)
# config -----------------------
outlier_file <- as.character(args[1])
outlier_regions <- read_tsv(outlier_file)

hap_to_perc <- 100 / (166 * 165)
iterations <- c("25/10 kb","10/5 kb","15/7.5 kb") %>% set_names(value = c("7", "8", "10"))
idx <- 10

import_map1 <- function(idx, filtmode = "direct"){
  read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_filterd.tsv")) %>% 
    mutate(ibd_total = (ibd2_cM_m1 + 0.5*ibd1_cM_m1) / (ibd0_cM_m1 + ibd1_cM_m1 + ibd2_cM_m1)) 
}

import_map2 <- function(idx, filtmode = "direct"){
  read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_filterd.tsv")) %>% 
    mutate(ibd_total = (ibd2_cM_m2 + 0.5*ibd1_cM_m2) / (ibd0_cM_m2 + ibd1_cM_m2 + ibd2_cM_m2)) 
}

import_bp <- function(idx, filtmode = "direct"){
  read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_filterd.tsv")) %>% 
    mutate(ibd_total = (ibd2_bp + 0.5*ibd1_bp) / (ibd0_bp + ibd1_bp + ibd2_bp)) 
}

plot_network <- function(idx, filt = 0, import_fun = import_map1, x = "cM_map1", filtmode = "direct"){
  clr2 <- GenomicOriginsScripts::clr[!(names(GenomicOriginsScripts::clr) %in% c("flo", "tor", "tab"))]
  clr2["uni"] <- rgb(.9,.9,.9)
  
  data <- import_fun(idx, filtmode = filtmode)
  
  set.seed(42)
  
  data %>% 
    as_tbl_graph() %E>%
    filter(ibd_total > filt) %N>%
    mutate(spec = str_sub(name,-6,-4),
           loc = str_sub(name,-3,-1))  %>% 
    ggraph( layout = 'fr', weights = ibd_total) +
    geom_edge_link(aes(alpha = ibd_total), color = rgb(.1,.1,.1), edge_width = .15) +
    geom_node_point(aes(fill = spec,
                        shape = loc, color = after_scale(clr_darken(fill,.3))), size = .7) +
    labs(y = glue::glue("Seq. Length: {iterations[as.character(idx)]}"),
         x = x) +
    scale_fill_manual("Species", values = GenomicOriginsScripts::clr[!(names(GenomicOriginsScripts::clr) %in% c("flo", "tor", "tab"))],
                      labels = GenomicOriginsScripts::sp_labs)+
    scale_edge_alpha_continuous(#limits = c(0,.1),
      range = c(0,1), guide = "none") +
    scale_shape_manual("Site", values = 21:23, labels = GenomicOriginsScripts::loc_names) +
    scale_x_continuous(position = "top") +
    guides(fill = guide_legend(title.position = "top",
                               nrow = 2, override.aes = list(shape = 21, size = 2.5)),
           shape = guide_legend(title.position = "top",
                                nrow = 2)) +
    coord_equal()  +
    theme(panel.background = element_blank(),
          axis.title.y = element_text(),
          axis.title.x = element_text())
}

p_done <- tibble(idx = rep(c(7,8,10), each = 3), 
       import_fun = rep(list(import_bp, import_map1, import_map2), 3),
       x = rep(c("bp", "cM_map1", "cM_map2"), 3)) %>% 
  pmap(plot_network) %>% 
  wrap_plots(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(text = element_text(size = plot_text_size),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.tag.position = c(0, 1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box = "horizontal", 
        legend.text.align = 0)

hypo_save(plot = p_done,
          filename = "figures/SFx1.png",
          width = .8*f_width,
          height = .8*f_width,
          dpi = 600,
          type = "cairo",
          bg = "transparent",
          comment = plot_comment)

p_filt <- tibble(idx = rep(c(7,8,10), each = 3), 
                 import_fun = rep(list(import_bp, import_map1, import_map2), 3),
                 x = rep(c("bp", "cM_map1", "cM_map2"), 3)) %>% 
  pmap(plot_network, filtmode = "bed95") %>% 
  wrap_plots(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(text = element_text(size = plot_text_size),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.tag.position = c(0, 1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box = "horizontal", 
        legend.text.align = 0)

hypo_save(plot = p_filt,
          filename = "figures/SFx1_filterd.png",
          width = .8*f_width,
          height = f_width,
          dpi = 600,
          type = "cairo",
          bg = "transparent",
          comment = plot_comment)
# ===========================
data_seg <- vroom::vroom(glue::glue("2_analysis/ibd/cM_converted/no_outgr_direct_10.conv_filterd.tsv"),
                         delim = "\t" ) %>%  #, col_types = "cccdddddddddddddddddddddd") %>% 
  mutate(ibd_hplo = str_remove(TYPE,"IBD") %>%
           as.integer())

y_lim <- c(0, 13)

p_bp <- data_seg %>%
  dplyr::select(seqnames = CHROM, start = bp_START, end = bp_END,G_SEG_START,GSTART,TYPE,ibd_hplo,PAIR) %>%
  arrange(G_SEG_START) %>% 
  dplyr::select(-G_SEG_START) %>% 
  as_granges() %>% 
  GenomicRanges::coverage(weight = "ibd_hplo") %>% 
  plyranges::as_ranges() %>% 
  as_tibble() %>% 
  dplyr::select(CHROM = seqnames, start, end, width, score) %>% 
  left_join(hypogen::hypo_chrom_start) %>% 
  mutate(gstart = GSTART + start, gend = GSTART + end) %>% 
  dplyr::select(CHROM, gstart, gend, score) %>% 
  pivot_longer(gstart:gend,values_to = "GPOS", names_to = "PART") %>% 
  ggplot() +
  geom_hypo_LG() +
  geom_vline(data = outlier_regions, aes(xintercept = gpos), color = rgb(1,0,0,.2), size = .3) +
  geom_step(aes(x = GPOS, y = score * hap_to_perc, group = CHROM), color = rgb(.3,.3,.3), size = .3) +
  geom_ribbon(aes(x = GPOS, ymin = 0, ymax = score * hap_to_perc, group = CHROM), fill = rgb(0,0,0,.5)) +
  scale_hypobg_manual(values = c("transparent",rgb(.9,.9,.9,.9),"red","blue") %>%
                        set_names(nm = c("even", "odd", "a","b")), guide = "none")+
  scale_x_hypo_LG() +
  labs(y = "IBD Score (bp)") +
  coord_cartesian(ylim = y_lim, expand = 0) +
  theme_hypo()

hypo_all_starts <- read_tsv(file = "ressources/hypo_all_starts.tsv") %>%
  mutate(grp = c("even", "odd")[row_number() %%2 + 1])

cM_digits <- 10^6

p_m1 <- data_seg %>%
  mutate(start = if_else(interpol_cM_START_m1 < interpol_cM_END_m1, interpol_cM_START_m1, interpol_cM_END_m1) * cM_digits,
         end = if_else(interpol_cM_START_m1 < interpol_cM_END_m1, interpol_cM_END_m1, interpol_cM_START_m1) * cM_digits) %>% 
  dplyr::select(seqnames = CHROM, start, end,G_SEG_START_cM_m1,GSTART_cM_m1,TYPE,ibd_hplo,PAIR) %>%
  arrange(G_SEG_START_cM_m1) %>% 
  dplyr::select(-G_SEG_START_cM_m1) %>% 
  as_granges() %>% 
  GenomicRanges::coverage(weight = "ibd_hplo") %>% 
  plyranges::as_ranges() %>% 
  as_tibble() %>% 
  dplyr::select(CHROM = seqnames, start, end, width, score) %>%
  mutate(start = start / cM_digits,
         end = end / cM_digits,
         width = width / cM_digits) %>% 
  left_join(hypo_all_starts) %>% 
  mutate(gstart = GSTART_cM_m1 + start, gend = GSTART_cM_m1 + end) %>% 
  dplyr::select(CHROM, gstart, gend, score) %>% 
  pivot_longer(gstart:gend,values_to = "GPOS", names_to = "PART") %>% 
  ggplot() +
  geom_rect(data = hypo_all_starts,
            aes(xmin = GSTART_cM_m1, xmax = GEND_cM_m1, ymin = -Inf, ymax = Inf, fill = grp)) +
  geom_step(aes(x = GPOS, y = score * hap_to_perc, group = CHROM), color = rgb(.3,.3,.3), size = .3) +
  geom_ribbon(aes(x = GPOS, ymin = 0, ymax = score * hap_to_perc, group = CHROM), fill = rgb(0,0,0,.5)) +
  coord_cartesian(ylim = y_lim, expand = 0) +
  scale_x_continuous(breaks = (hypo_all_starts$GSTART_cM_m1 + hypo_all_starts$GEND_cM_m1)/2,
                     labels = hypo_all_starts$CHROM %>%
                       str_remove("LG"),
                     position = "top",
                     limits = c(0, max(hypo_all_starts$GEND_cM_m1)),expand = c(0, 0)) +
  scale_fill_manual(values = c(odd = rgb(.6,.6,.6,.3), even = "transparent"), guide = FALSE) +
  labs(y = "IBD Score (cM1)") +
  theme_hypo()


p_m2 <- data_seg %>%
  mutate(start = if_else(interpol_cM_START_m2 < interpol_cM_END_m2, interpol_cM_START_m2, interpol_cM_END_m2) * cM_digits,
         end = if_else(interpol_cM_START_m2 < interpol_cM_END_m2, interpol_cM_END_m2, interpol_cM_START_m2) * cM_digits) %>% 
  dplyr::select(seqnames = CHROM, start, end,G_SEG_START_cM_m2,GSTART_cM_m2,TYPE,ibd_hplo,PAIR) %>%
  arrange(G_SEG_START_cM_m2) %>% 
  dplyr::select(-G_SEG_START_cM_m2) %>% 
  as_granges() %>% 
  GenomicRanges::coverage(weight = "ibd_hplo") %>% 
  plyranges::as_ranges() %>% 
  as_tibble() %>% 
  dplyr::select(CHROM = seqnames, start, end, width, score) %>%
  mutate(start = start / cM_digits,
         end = end / cM_digits,
         width = width / cM_digits) %>% 
  left_join(hypo_all_starts) %>% 
  mutate(gstart = GSTART_cM_m2 + start, gend = GSTART_cM_m2 + end) %>% 
  dplyr::select(CHROM, gstart, gend, score) %>% 
  pivot_longer(gstart:gend,values_to = "GPOS", names_to = "PART") %>% 
  ggplot() +
  geom_rect(data = hypo_all_starts,
            aes(xmin = GSTART_cM_m2, xmax = GEND_cM_m2, ymin = -Inf, ymax = Inf, fill = grp)) +
  geom_step(aes(x = GPOS, y = score * hap_to_perc, group = CHROM), color = rgb(.3,.3,.3), size = .3) +
  geom_ribbon(aes(x = GPOS, ymin = 0, ymax = score * hap_to_perc, group = CHROM), fill = rgb(0,0,0,.5)) +
  coord_cartesian(ylim = y_lim, expand = 0) +
  scale_x_continuous(breaks = (hypo_all_starts$GSTART_cM_m2 + hypo_all_starts$GEND_cM_m2)/2,
                     labels = hypo_all_starts$CHROM %>%
                       str_remove("LG"),
                     position = "top",
                     limits = c(0, max(hypo_all_starts$GEND_cM_m2)), expand = c(0, 0)) +
  scale_fill_manual(values = c(odd = rgb(.6,.6,.6,.3), even = "transparent"), guide = FALSE) +
  labs(y = "IBD Score (cM2)") +
  theme_hypo()

p_gw <- p_bp /
  p_m1 /
  p_m2 & 
  theme(text = element_text(size = plot_text_size))

hypo_save(plot = p_gw,
          filename = "figures/SFx1_gw_filtered.pdf",
          width = f_width,
          height = .6*f_width,
          device = cairo_pdf,
          # dpi = 600,
          # type = "cairo",
          bg = "transparent",
          comment = plot_comment)
