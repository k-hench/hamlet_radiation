#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFX2.R 2_analysis/ccf/
# ===============================================================
# This script produces Suppl. Figure X of the study "Ancestral variation,
# hybridization and modularity fuel a marine radiation"
# by Hench, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "2_analysis/ccf/" )
# script_name <- "R/fig/plot_SFX.R"
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(GenomicOriginsScripts)
library(ggtext)
library(scales)
library(grid)

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
ccf_dir <- as.character(args[1])

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
 
read_ccf <- function(path = ccf_dir, file ){
  ind1 <- file %>% str_remove("^ccf.LG[0-9]{2}.") %>% str_split(pattern = "\\.") %>% .[[1]] %>% .[[1]]
  ind2 <- file %>% str_remove("^ccf.LG[0-9]{2}.") %>% str_split(pattern = "\\.") %>% .[[1]] %>% .[[2]]
  vroom::vroom(str_c(path,file)) %>%
    set_names(nm = c("CHROM", "POS", "PostMedian",
                     str_c(ind1, "_1-",ind1, "_2"),
                     str_c(ind1, "_1-",ind2, "_1"),
                     str_c(ind1, "_1-",ind2, "_2"),
                     str_c(ind1, "_2-",ind2, "_1"),
                     str_c(ind1, "_2-",ind2, "_2"),
                     str_c(ind2, "_2-",ind2, "_2"))) %>%
    pivot_longer(names_to = "pair", 
                 cols = -(CHROM:PostMedian)) %>%
    mutate(prep = pair) %>%
    separate(prep, into = c("ind1","ind2"), sep = "-") %>%
    separate(ind1, into = c("ind1","hap1"), sep = "_",convert = TRUE)%>%
    separate(ind2, into = c("ind2","hap2"), sep = "_",convert = TRUE) %>%
    mutate(spec1 = ind1 %>% str_sub(-6,-4),
           spec2 = ind2 %>% str_sub(-6,-4),
           value = if_else(hap1 == 1, value, -value)) %>%
    arrange(pair, PostMedian) %>%
    group_by(pair) %>%
    mutate(check = lag(value, default = 0) == value & lead(value, default = 0) == value) %>%
    ungroup() %>%
    filter(!check)
}

dummy <- tibble(hap1 = 1:2,
                start = c(0,0),
                end = c(1,-1))

xrange <- c(.1,10^6)
color <- rgb(1, 0.5, 0.16)

base_length <- 8
base_lwd <- .15
base_line_clr <- "black"

plot_ccf <- function(lg, target, ccf_dir_in = ccf_dir){
  
  ccf_files <- dir(ccf_dir_in, pattern = str_c("ccf.", lg, ".",target))
  print(ccf_files)
  
  if(length(ccf_files > 0)){
  ccf_data <- ccf_files %>% 
    purrr::map_dfr(read_ccf, path = ccf_dir) 
  
  p <- ccf_data  %>%
    ggplot(aes(x = PostMedian, y = value, group = pair, color = spec2)) +
    geom_segment(inherit.aes = FALSE,
                 data = dummy, aes(x = 10, xend = 10, y = start, yend = end),
                 color =rgb(1,1,1,0))+
    theme_gradient_bg_x(xlim_in = xrange, #range(economics$unemploy),
                        colors = c( rgb(.95,.95,.95,0),
                                    rgb(.1,.1,.1,.15)),
                        fun = log10) +
    geom_step()+
    facet_grid(hap1 ~ ., scales = "free_y")+
    labs(y = "Fraction of genome shared", title = str_c(lg," / ", target))+
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_color_manual(values = clr2)+
    theme_minimal()+
    theme(text = element_text(family = "Futura Lt BT"),
          axis.title.y = element_markdown(),
          legend.position = "bottom",
          plot.subtitle = element_markdown(),
          axis.line = element_line(colour = base_line_clr,
                                   size = base_lwd), 
          strip.background = element_rect(fill = rgb(.95,.95,.95),colour = base_line_clr,size = base_lwd),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  
  ggsave(plot = p,
         filename = str_c("~/Desktop/tests/figdump/SFX2_",lg,"_",target,".pdf"),
         width = 7,height = 4, device = cairo_pdf, bg = "transparent")}
  else {print( str_c("skip ", lg, " ",target, "\\n"))}
}


cross_df(list(lg = c("LG04", "LG06","LG12", "LG15", "LG17"),
              target  = c("18152puebel", "18161puebel", "18172puebel", "18174puebel", "18180puebel"))) %>%
  purrr::pmap(plot_ccf)

