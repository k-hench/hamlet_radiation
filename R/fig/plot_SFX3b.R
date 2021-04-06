#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFX3.R
# ===============================================================
# This script produces Suppl. Figure X of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "2_analysis/summaries/fst_outliers_998.tsv", "2_analysis/geva/", "2_analysis/GxP/bySNP/" , "2_analysis/GxP/50000/")
# script_name <- "R/fig/plot_SFX3b.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(hypogen)
library(hypoimg)
library(ggtext)
library(ggpointdensity)
library(scales)
library(grid)
library(geomfactory)
library(prismatic)
library(patchwork)

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
outlier_file <- as.character(args[1])
geva_path <- as.character(args[2])
gxp_path <- as.character(args[3])
gxp_win_path <- as.character(args[4])

factory_geom_line('gxp')
factory_geom_point('gxp')
#factory_geom_point('snpgxp')
outlier_data <- read_tsv(outlier_file)

window_buffer <- 2.5*10^5

import_geva_data <- function(chrom, start, end, ...){
  vroom::vroom(file = str_c(geva_path, chrom,".sites.txt.gz"), delim = " ") %>%
    left_join(vroom::vroom(str_c(geva_path, chrom,".marker.txt.gz"), delim = " ")) %>%
    arrange(Position) %>%
    filter(between(Position,start-window_buffer,end+window_buffer)) %>%
    select(Chromosome,Position,MarkerID, Clock, Filtered:PostMedian) %>%
    mutate(CHROM = str_c("LG", str_pad(Chromosome,width = 2,pad = "0")))
}

gxp_importer <- function(trait, chrom, start, end ){
  vroom::vroom(str_c(gxp_path,trait,".lm.GxP.txt.gz"), delim = "\t") %>%
    filter(CHROM == chrom, between(POS, start-window_buffer, end+window_buffer) ) %>%
    select(CHROM, POS, p_wald) %>%
    set_names(nm = c("CHROM", "Position", trait))
}

import_gxp_data <- function(chrom, start, end, ...){
  bar_data <- gxp_importer("Bars", chrom, start, end)
  peduncle_data <- gxp_importer("Peduncle", chrom, start, end)
  snout_data <-gxp_importer("Snout", chrom, start, end)
  
  bar_data %>%
    left_join(peduncle_data, by = c(CHROM = "CHROM", Position = "Position")) %>%
    left_join(snout_data, by = c(CHROM = "CHROM", Position = "Position"))
}

get_gxp_and_geva <- function(gid, chrom, start, end, ...){
  import_geva_data(chrom, start, end) %>%
    left_join(import_gxp_data(chrom, start, end)) %>%
    mutate(gid = gid)
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

grid_piece_x <- function(g, x1, x2, ...){
  annotation_custom(grob = g, xmin = x1, xmax = x2, ...)
}

theme_gradient_bg_x <- function(xlim_in, colors = c("#7B0664", "#E32219"), fun = identity, step = 1,...){
  g <- rasterGrob(cbind(colors[[1]], colors[[2]]), width = unit(1, "npc"),
                  height = unit(1, "npc"), 
                  interpolate = TRUE)
  
  xstart <- floor(fun(min(xlim_in)))
  xend <- ceiling(fun(max(xlim_in)))
  
  breaks <- seq(xstart,xend, by = step)
  
  purrr::map2(.x = breaks[1:(length(breaks)-1)],
              .y = breaks[2:length(breaks)],
              .f = grid_piece_x, 
              g = g, ...)
}

trimm_to_outlier <- function(data, gid, chrom, start, end, ...){
  data %>%
    filter(CHROM == chrom, between(MID_POS, start-window_buffer, end+window_buffer)) %>%
    mutate(gid = gid)
}


outlier_gxp_data <- function(trait, outlier_data){
  
  outlier_data %>%
    pmap_dfr( trimm_to_outlier, data = vroom::vroom(str_c("2_analysis/GxP/50000/",trait,".lm.50k.5k.txt.gz"))) %>%
    mutate(trait = trait)
}


get_anno <- function (gid, chrom, start, end, anno_rown = 5, ...) {
  df_list <- hypo_annotation_get(searchLG = chrom, xrange = c(start - window_buffer, end + window_buffer), anno_rown = anno_rown)
  df_list[[1]] <- df_list[[1]] %>% mutate(gid = gid)
  df_list[[2]] <- df_list[[2]] %>% mutate(gid = gid)
  
  tibble(exnons = list(df_list[[2]]), genes = list(df_list[[1]] ))
}

gxp_win_data <- list("Bars", "Snout", "Peduncle") %>%
  map_dfr(.f = outlier_gxp_data, outlier_data  = outlier_data[c(2,13,14),]) %>%
  mutate(POS = MID_POS, window = "GxP  (-log<sub>10</sub>  *p*-value)")

annos <- outlier_data[c(2,13,14),] %>%
  pmap_dfr(get_anno)

exons <- bind_rows(annos$exnons) %>% mutate(window = "genes")
genes <-  bind_rows(annos$genes) %>% mutate(window = "genes")

data <- outlier_data[c(2,13,14),] %>%
  pmap_dfr(get_gxp_and_geva)

xrange <- c(100,10^6)
color <- rgb(1, 0.5, 0.16)

base_length <- 8
base_lwd <- .15
base_line_clr <- "black"
base_fnt <-  "Futura Lt BT"

plot_data <- data %>%
  pivot_longer(names_to = "trait",
               values_to = "p_wald",
               cols = Bars:Snout) %>%
  filter(Clock == "J", Filtered == 1) %>%
  mutate(POS = Position, CHROM = Chromosome, window = "PosteriorMedian (log<sub>10</sub>)")

splitage <- tibble(intercept = log10(5000), window = "PosteriorMedian (log<sub>10</sub>)")

gxp_clr <- c(Bars = "#79009f", Snout = "#E48A00", Peduncle = "#5B9E2D") %>%
  darken(factor = .95) %>%
  set_names(., nm = c("Bars", "Snout", "Peduncle"))

p <- plot_data %>%
  ggplot(aes(x = POS)) +
  facet_grid(window~gid, scales = "free")+
  ggplot2::geom_segment(data = (genes %>% filter(strand %in% c("+", "-"))), 
                        aes(x = ps, xend = pe, y = yl, yend = yl, group = Parent), 
                        lwd = base_lwd, arrow = arrow(length = unit(2, "pt"), type = "closed"),
                        color = base_line_clr) +
  ggplot2::geom_segment(data = (genes %>%  filter(!strand %in% c("+", "-"))), 
                        aes(x = ps, xend = pe,  y = yl, yend = yl, group = Parent),
                        lwd = base_lwd, color = base_line_clr) + 
  ggplot2::geom_text(data = genes, size = 3, 
                     aes(x = labelx, label = gsub("hpv1g000000",  ".", label), y = yl - 0.5),
                     #family = base_fnt, 
                     fontface = "italic")+
  geom_line_gxp(data = gxp_win_data, aes( y = AVG_p_wald, gxp = trait))+
  geom_pointdensity(aes(y = log10(PostMedian), color = ..density..), size = .85)+
  geom_point_gxp( data = plot_data %>%
                    filter(trait == "Bars" & gid == "LG12_3" |
                             trait == "Peduncle" & gid == "LG12_4" |
                             trait == "Snout" & gid == "LG04_1") %>%
                    mutate(window = "GxP by SNP (-log<sub>10</sub>  *p*-value)"),
                  aes(y = -log10(p_wald), gxp_c =trait), 
                  size = .5, alpha = .2)+
  geom_hline(data = splitage, aes(yintercept = intercept), 
             linetype = 3, color = rgb(1,1,1,.9)) +
  scale_color_viridis_c(option = "B",
                        labels = scales::trans_format("log10", scales::math_format(10^.x)), breaks = c(10^-6, 10^-5.5))+
  scale_gxp_manual(values = gxp_clr) +
  scale_gxp_c_manual(values = gxp_clr)+
  # scale_snpgxp_c_continuous()+#low = "black",high = "darkred")+
  scale_alpha_continuous(guide = FALSE) +
  guides(color = guide_colorbar(barwidth = unit(150, "pt"),
                                barheight = unit(7, "pt"))#,
         #snpgxp_c = guide_colourbar_snpgxp(barwidth = unit(150, "pt"),
         #                       barheight = unit(7, "pt"))
  )+
  theme_minimal()+
  theme(#text = element_text(family = base_fnt),
    axis.title.y = element_markdown(),
    legend.position = "bottom",
    plot.subtitle = element_markdown(),
    axis.line = element_line(colour = base_line_clr,
                             size = base_lwd), 
    strip.background = element_rect(fill = rgb(.95,.95,.95),colour = base_line_clr,size = base_lwd),
    panel.grid.minor = element_blank(),
    strip.text.y  = element_markdown())

subplot()

scl <- 1
hypo_save(plot = p,
          filename = "figures/SFX3.pdf",
          width = 12 * scl,height = 8 * scl,
          comment = plot_comment,
          device = cairo_pdf,
          bg = "transparent")

# ------------------------------------------------------------------------------

theme_subplot <- function(...){
  theme_minimal()+
    theme(axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
}

plot_panel2_anno <- function(gid_in = "LG12_3", title = ""){
  genes_of_interest <-  c('arl3','kif16b','cdx1','hmcn2',
                                        'sox10','smarca4',
                                        'rorb',
                                        'alox12b','egr1',
                                        'ube4b','casz1',
                                        'hoxc8a','hoxc9','hoxc10a',
                                        'hoxc13a','rarga','rarg',
                                        'snai1','fam83d','mafb','sws2abeta','sws2aalpha','sws2b','lws','grm8')
  outlier_data %>%
    filter(gid == gid_in) %>%
    ggplot() +
    coord_cartesian(xlim = c(outlier_data$start[[which(outlier_data$gid == gid_in)]] - window_buffer,
                             outlier_data$end[[which(outlier_data$gid == gid_in)]] + window_buffer),
                    expand = .01,
                    ylim = c(1,6.5))+
    geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), color = rgb(0,0,0,.2), fill = rgb(0,0,0,.1)) +
    ggplot2::geom_segment(data = (genes %>%
                                    filter(gid == gid_in) %>% 
                                    filter(strand %in% c("+", "-"))), 
                          aes(x = ps, xend = pe, y = yl, yend = yl, group = Parent), 
                          lwd = base_lwd, arrow = arrow(length = unit(2, "pt"), type = "closed"),
                          color = base_line_clr) +
    ggplot2::geom_segment(data = (genes %>%
                                    filter(gid == gid_in) %>%
                                    filter(!strand %in% c("+", "-"))), 
                          aes(x = ps, xend = pe,  y = yl, yend = yl, group = Parent),
                          lwd = base_lwd, color = base_line_clr) + 
    ggplot2::geom_text(data = genes %>%
                         filter(gid == gid_in,
                                label %in% genes_of_interest),
                       size = 3, 
                       aes(x = labelx, label = gsub("hpv1g000000",  ".", label), y = yl - 0.5),
                       #family = base_fnt, 
                       fontface = "italic")+
    scale_y_continuous("")+
    labs(title = title)+
    theme_subplot()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(hjust = .5))}

plot_panel2_gxp <- function(gid_in = "LG12_3"){
    ggplot(mapping = aes(x = POS)) +
    coord_cartesian(xlim = c(outlier_data$start[[which(outlier_data$gid == gid_in)]] - window_buffer,
                             outlier_data$end[[which(outlier_data$gid == gid_in)]] + window_buffer),
                    expand = .01)+
    geom_hypo_grob(inherit.aes = FALSE,
                   data = tibble(grob = list(annotation_grobs[[gid_in]])), 
                   aes(x = .82, y = .75, grob = grob), width = .42, height = .42, angle = 0)+
    geom_rect(data = outlier_data %>%
                filter(gid == gid_in),
              inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              color = rgb(0,0,0,.2), fill = rgb(0,0,0,.1)) +
    geom_line(data = gxp_win_data%>%
                filter(gid == gid_in), aes( y = AVG_p_wald, color = trait))+
    scale_color_manual(values = gxp_clr)+
    scale_y_continuous("GxP<sub>50 kb</sub> (-log<sub>10</sub>  *p*)")+
    theme_subplot()
}

plot_panel2_gxp_point <- function(gid_in = "LG12_3"){
    ggplot(mapping = aes(x = POS)) +
    coord_cartesian(xlim = c(outlier_data$start[[which(outlier_data$gid == gid_in)]] - window_buffer,
                             outlier_data$end[[which(outlier_data$gid == gid_in)]] + window_buffer),
                    expand = .01)+
    geom_rect(data = outlier_data %>%
                filter(gid == gid_in),
              inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              color = rgb(0,0,0,.2), fill = rgb(0,0,0,.1)) +
    geom_point( data = plot_data %>%
                  filter(gid == gid_in)%>%
                      filter(trait == "Bars" & gid == "LG12_3" |
                               trait == "Peduncle" & gid == "LG12_4" |
                               trait == "Snout" & gid == "LG04_1") %>%
                      mutate(window = "GxP by SNP (-log<sub>10</sub>  *p*-value)"),
                    aes(y = -log10(p_wald), color =trait), 
                    size = .5, alpha = .2)+
    scale_color_manual(values = gxp_clr, guide = FALSE)+
    scale_y_continuous("GxP<sub>SNP</sub> (-log<sub>10</sub>  *p*)")+
    theme_subplot()
}

plot_panel2_aa <- function(gid_in = "LG12_3"){
    ggplot(mapping = aes(x = POS)) +
    coord_cartesian(xlim = c(outlier_data$start[[which(outlier_data$gid == gid_in)]] - window_buffer,
                             outlier_data$end[[which(outlier_data$gid == gid_in)]] + window_buffer),
                    expand = .01,ylim = c(10^1.6,10^5.6))+
    geom_rect(data = outlier_data %>%
                filter(gid == gid_in),
              inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 10^1.6, ymax = 10^5.6),
              color = rgb(0,0,0,.2), fill = rgb(0,0,0,.1)) +
    geom_pointdensity(data = plot_data%>%
                        filter(gid == gid_in),aes(y = PostMedian, color = ..density..), size = .85)+
    # scale_color_viridis_c(option = "D", 
    scico::scale_color_scico(palette = "berlin",
                             limits = c(10^-10,10^-5.3), 
                             labels = scales::trans_format("log10", scales::math_format(10^.x)),
                             breaks = c(10^-10, 10^-5.5))+
    scale_alpha_continuous(guide = FALSE) +
    scale_y_log10("Allele Age (Gen.)",
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    guides(color = guide_colorbar(title = "Density",
                                  title.position = "top",
                                  barwidth = unit(150, "pt"),
                                  barheight = unit(7, "pt"),
                                  direction = "horizontal")
    )+
    theme_minimal()
}

plot_curtain2 <- function(gid_in, title, theme_curtain = theme){
            plot_grid(plot_panel2_anno(gid_in, title) + theme_curtain(), 
            plot_panel2_gxp(gid_in) + theme_curtain(), 
            plot_panel2_gxp_point(gid_in) + theme_curtain(),
            plot_panel2_aa(gid_in) + theme_curtain(),
            ncol = 1,align = "v",label_fontface = "plain")
}

annotation_grobs <- tibble(svg = hypo_trait_img$grob_circle[hypo_trait_img$trait %in% c( 'Snout', 'Bars', 'Peduncle')],
                          layer = c(4,3,7),
                          color = gxp_clr[c(1,3,2)]) %>%
  pmap(hypo_recolor_svg) %>%
  set_names(nm = c( "LG12_3","LG12_4","LG04_1"))
annotation_grobs$LG12_3 <- hypo_recolor_svg(annotation_grobs$LG12_3,layer = 7, color = gxp_clr[[1]] %>% clr_desaturate %>% clr_lighten(.25))


theme_inner_curtain <- function(...){theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), ...)}
theme_outer_curtain <- function(...){theme(legend.position = "none", axis.title.y = element_markdown(), axis.title.x = element_blank(), ...)}

p_prep <- (plot_curtain2("LG04_1", title = "LG04_1 (A)", theme_curtain = theme_outer_curtain ) | 
  plot_curtain2("LG12_3", title = "LG12_3 (B)", theme_curtain = theme_inner_curtain) |
plot_curtain2("LG12_4", title = "LG12_4 (C)", theme_curtain = theme_inner_curtain)) +
  plot_annotation(tag_levels = "a")

leg_gxp <- (plot_panel2_gxp() + guides(color = guide_legend(title = "Trait",title.position = "top")) + theme(legend.position = "bottom")) %>% get_legend()
leg_aa <- (plot_panel2_aa() + theme(legend.position = "bottom")) %>% get_legend()

p <- plot_grid(p_prep, 
          plot_grid(leg_gxp, leg_aa, nrow = 1),
          ncol = 1,rel_heights = c(1,.1))

scl <- 1
hypo_save(plot = p,
          filename = "figures/SFX3b.pdf",
          width = 12 * scl,height = 8 * scl,
          comment = plot_comment,
          device = cairo_pdf,
          bg = "transparent")
