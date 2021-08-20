#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF13.R \
#     2_analysis/summaries/fst_outliers_998.tsv
# ===============================================================
# This script produces Suppl. Figure 13 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c("2_analysis/summaries/fst_outliers_998.tsv")
# script_name <- "R/fig/plot_SF13.R"
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
iterations <- c(str_c("2/5*10^",6:3," BP"),"25/10 kb","10/5 kb","7-5/3 kb","15/7.5 kb")
itteration_names <-   c(str_c("10-",6:3),"7","8","9","10")

plot_network <- function(n_zeros, filtmode = "direct", filt = 0){
  clr2 <- GenomicOriginsScripts::clr[!(names(GenomicOriginsScripts::clr) %in% c("flo", "tor", "tab"))]
  clr2["uni"] <- rgb(.9,.9,.9)
  
  data <- read_tsv(glue::glue("2_analysis/ibd/no_outgr_{filtmode}_{itteration_names[n_zeros]}.ibd.tsv")) %>% 
    mutate(ibd_total = (IBD2 + 0.5*IBD1) / (IBD0 + IBD1 + IBD2)) 
  
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
    labs(y = glue::glue("Seq. Length: {iterations[n_zeros]}")) +
    scale_fill_manual("Species", values = GenomicOriginsScripts::clr[!(names(GenomicOriginsScripts::clr) %in% c("flo", "tor", "tab"))],
                      labels = GenomicOriginsScripts::sp_labs)+
    scale_edge_alpha_continuous(#limits = c(0,.1),
                                range = c(0,1), guide = "none") +
    scale_shape_manual("Site", values = 21:23, labels = GenomicOriginsScripts::loc_names) +
    guides(fill = guide_legend(nrow = 2, override.aes = list(shape = 21, size = 2.5)),
           shape = guide_legend(nrow = 2)) +
    coord_equal()  +
    theme(panel.background = element_blank(),
          axis.title.y = element_text())
}

plot_ibd_gw <- function(n_zeros, y_lim = c(0, 4), filtmode = "direct"){
  data_seg <- vroom::vroom(glue::glue("2_analysis/ibd/no_outgr_{filtmode}_{itteration_names[n_zeros]}.segments.tsv"), delim = "\t", col_types = "cccciidcdci") %>% 
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
    labs(y = "IBD Score") +
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

plts <- c(c(5,8,6) %>% map(plot_network),
          c(5,8,6) %>% map(plot_network, filtmode = "bed95") %>% map(.f = function(p){ p + theme(axis.title.y = element_blank())}),
          c(5,8,6) %>% map2(.y = list(c(0,26), # .55),
                                      c(0,26), # 4),
                                      c(0,26)),
                            plot_ibd_gw,  filtmode = "direct"
          ))

p_done <- (plts[[1]] + plts[[4]] + plts[[7]] + 
             plts[[2]] + plts[[5]] + plts[[8]] + geom_hline(yintercept = perc_score * hap_to_perc, color = "#11C269", size = .3, alpha = .7) + 
             plts[[3]] + plts[[6]] + plts[[9]] + 
             plot_layout(ncol = 3,
                         widths = c(.3,.3, 1)))/ 
  guide_area() +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(1, .07),
              guides = "collect") &
  theme(text = element_text(size = plot_text_size),
        plot.tag.position = c(0, 1),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box = "horizontal", 
        legend.text.align = 0)

hypo_save(plot = p_done,
       filename = "figures/SF13.png",
       width = f_width,
       height = .7 * f_width,
       dpi = 600,
       type = "cairo",
       bg = "transparent",
       comment = plot_comment)

system("convert figures/SF13.png figures/SF13.pdf")
system("rm figures/SF13.png")
create_metadata <- str_c("exiftool -overwrite_original -Description=\"", plot_comment, "\" figures/SF13.pdf")
system(create_metadata)