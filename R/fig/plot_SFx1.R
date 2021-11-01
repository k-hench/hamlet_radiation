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

import_map1 <- function(idx, filtmode = "bed95"){
  read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_filterd.tsv")) %>% 
    mutate(ibd_total = (ibd2_cM_m1 + 0.5*ibd1_cM_m1) / (ibd0_cM_m1 + ibd1_cM_m1 + ibd2_cM_m1)) 
}

import_map2 <- function(idx, filtmode = "bed95"){
  read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_filterd.tsv")) %>% 
    mutate(ibd_total = (ibd2_cM_m2 + 0.5*ibd1_cM_m2) / (ibd0_cM_m2 + ibd1_cM_m2 + ibd2_cM_m2)) 
}

import_bp <- function(idx, filtmode = "bed95"){
  # read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_filterd.tsv")) %>% 
    read_tsv(glue::glue("2_analysis/ibd/cM_converted/no_outgr_{filtmode}_{idx}.conv_summary.tsv")) %>% 
    mutate(ibd_total = (ibd2_bp + 0.5*ibd1_bp) / (ibd0_bp + ibd1_bp + ibd2_bp)) 
}

import_truffle <- function(idx, filtmode = "direct"){
  itteration_names <- c(str_c("10-",6:3),"7","8","9","10")
  read_tsv(glue::glue("2_analysis/ibd/no_outgr_{filtmode}_{itteration_names[idx]}.ibd.tsv")) %>% 
    mutate(ibd_total = (IBD2 + 0.5*IBD1) / (IBD0 + IBD1 + IBD2)) 
}

iterations <- c(str_c("2/5*10^",6:3," BP"),"25/10 kb","10/5 kb","7-5/3 kb","15/7.5 kb")

plot_network <- function(idx, filt = 0, import_fun = import_map1,
                         x = "cM_map1", filtmode = "direct", 
                         x_ax = TRUE, y_ax = TRUE, ...){
  clr2 <- GenomicOriginsScripts::clr[!(names(GenomicOriginsScripts::clr) %in% c("flo", "tor", "tab"))]
  clr2["uni"] <- rgb(.9,.9,.9)
  
  data <- import_fun(idx, filtmode = filtmode)
  
  set.seed(42)
  
  p <- data %>% 
    as_tbl_graph() %E>%
    filter(ibd_total > filt) %N>%
    mutate(spec = str_sub(name,-6,-4),
           loc = str_sub(name,-3,-1))  %>% 
    ggraph( layout = 'fr', weights = ibd_total) +
    geom_edge_link(aes(alpha = ibd_total), color = rgb(.1,.1,.1), edge_width = .15) +
    geom_node_point(aes(fill = spec,
                        shape = loc, color = after_scale(clr_darken(fill,.3))), size = .7) +
    labs(y = glue::glue("Seq. Length: {iterations[idx]}"),
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
                                nrow = 2, override.aes = list(size = 2.5))) +
    coord_equal()  +
    theme(panel.background = element_blank(),
          axis.title.y = element_text(),
          axis.title.x = element_text())
  
  if(!x_ax){ p <- p + theme(axis.title.x = element_blank())}
  if(!y_ax){ p <- p + theme(axis.title.y = element_blank())}
  p
}

plts_cM <- tibble(idx = rep(c(7, 10, 8), each = 3), 
       import_fun = rep(list(import_bp, import_map1, import_map2), 3),
       x = rep(c("bp_cM_filt.", "cM_map1", "cM_map2"), 3),
       x_ax = rep(c(TRUE, FALSE), c(3, 6)),
       y_ax = rep(FALSE, 9),
       filtmode = "bed95") %>% 
  bind_rows(tibble(idx = rep(c(5, 8, 6), 2), 
                   import_fun = rep(list(import_truffle), 6),
                   x = rep(c("truffle", "bed95"), each = 3),
                   x_ax = rep(rep(c(TRUE, FALSE), 1:2), 2),
                   y_ax = rep(c(TRUE, FALSE), each = 3),
                   filtmode = rep(c("direct", "bed95"), each = 3)) ) %>% 
  left_join(tibble(x = c("truffle", "bed95", "bp_cM_filt.", "cM_map1", "cM_map2"),
            plot_order = seq_along(x))) %>% 
  arrange(plot_order) %>% 
  pmap(plot_network) 

p_done <- plts_cM %>% 
  wrap_plots(nrow = 3,
             byrow = FALSE,
             guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(text = element_text(size = plot_text_size),
        plot.tag.position = c(0, 1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box = "horizontal", 
        legend.text.align = 0,
        plot.subtitle = element_text())

hypo_save(plot = p_done,
          filename = "figures/SFxx1.png",
          width = f_width,
          height = .75*f_width,
          dpi = 600,
          type = "cairo",
          bg = "transparent",
          comment = plot_comment)

system("convert figures/SFxx1.png figures/SFxx1.pdf")
system("rm figures/SFxx1.png")
create_metadata <- str_c("exiftool -overwrite_original -Description=\"", plot_comment, "\" figures/SFxx1.pdf")
system(create_metadata)

# -------------------------------------------------
plot_ibd_gw <- function(n_zeros, y_lim = c(0, 4), filtmode = "direct", y_lab = ""){
  data_seg <- vroom::vroom(glue::glue("2_analysis/ibd/no_outgr_{filtmode}_{itteration_names[n_zeros]}.segments.tsv"),
                           delim = "\t", col_types = "cccciidcdci") %>% 
    left_join(hypogen::hypo_chrom_start) %>% 
    mutate(start = POS * 10^6,
           end = start + (LENGTH * 10^6),
           gstart = GSTART + start,
           gend = start + (LENGTH * 10^6),
           ibd_hplo = str_remove(TYPE,"IBD") %>%
             as.integer())
  data_seg %>%
    dplyr::select(seqnames = CHROM,start,end,gstart,GSTART,TYPE,ibd_hplo,ID1,ID2) %>%
    arrange(gstart) %>% 
    dplyr::select(-gstart) %>% 
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
    labs(y = str_c("IBD Score (" , iterations[n_zeros], ")")) +
    coord_cartesian(ylim = y_lim, expand = 0) +
    theme_hypo()
}

data_seg <- vroom::vroom(glue::glue("2_analysis/ibd/no_outgr_direct_10.segments.tsv"),
                         delim = "\t", col_types = "cccciidcdci") %>%
  left_join(hypogen::hypo_chrom_start) %>%
  mutate(start = POS * 10^6,
         end = start + (LENGTH * 10^6),
         gstart = GSTART + start,
         gend = start + (LENGTH * 10^6),
         ibd_hplo = str_remove(TYPE,"IBD") %>%
           as.integer())

data_compact <- data_seg %>%
  dplyr::select(seqnames = CHROM,start,end,gstart,GSTART,TYPE,ibd_hplo,ID1,ID2) %>%
  arrange(gstart) %>%
  dplyr::select(-gstart) %>%
  as_granges() %>%
  GenomicRanges::coverage(weight = "ibd_hplo") %>%
  plyranges::as_ranges() %>%
  as_tibble() %>%
  dplyr::select(CHROM = seqnames, start, end, width, score) %>%
  left_join(hypogen::hypo_chrom_start) %>%
  mutate(gstart = GSTART + start, gend = GSTART + end) %>%
  dplyr::select(CHROM, gstart, gend, score)

total_cov_lenght <- data_compact %>%
  mutate(length = gend-gstart) %>% .$length %>% sum()

data_sorted <- data_compact %>%
  mutate(length = gend-gstart) %>%
  group_by(score) %>%
  summarise(length = sum(length)) %>%
  ungroup() %>%
  mutate(length = if_else(score == 0,
                          length + hypo_karyotype$GEND[24]-total_cov_lenght, # attach uncovered chrom ends
                          length),
         gend = cumsum(length),
         gstart = lag(gend,default = 0))

perc_cutoff <- .95

perc_score <- data_sorted %>%
  filter(gstart < hypo_karyotype$GEND[24] * perc_cutoff,
         gend > hypo_karyotype$GEND[24] * perc_cutoff) %>%
  .$score

# tibble(n_zeros = c(5,8,6),
#        y_lim = list(c(0,26), # .55),
#                     c(0,26), # 4),
#                     c(0,26)),
#        ylab = 

plts <- c(5,8,6) %>%
  map2(.y = list(c(0,26), # .55),
                 c(0,26), # 4),
                 c(0,26)),
       plot_ibd_gw, filtmode = "direct")

p_done <- (plts[[1]] +
             plts[[2]] + geom_hline(yintercept = perc_score * hap_to_perc, color = "#11C269", size = .3, alpha = .7) + 
             plts[[3]] +
             plot_layout(ncol = 1)) +
  plot_annotation(tag_levels = "a") &
  theme(text = element_text(size = plot_text_size),
        plot.tag.position = c(0, 1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box = "horizontal", 
        legend.text.align = 0)

hypo_save(plot = p_done,
          filename = "figures/SFxx2.png",
          width = f_width,
          height = .55 * f_width,
          dpi = 600,
          type = "cairo",
          bg = "transparent",
          comment = plot_comment)

system("convert figures/SFxx2.png figures/SFxx2.pdf")
system("rm figures/SFxx2.png")
create_metadata <- str_c("exiftool -overwrite_original -Description=\"", plot_comment, "\" figures/SFxx2.pdf")
system(create_metadata)
