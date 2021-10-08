#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/recombination_scaling.R \
#     ressources/recombination/ \
#     MAP1cm.txt MAP1bp.txt \
#     MAP2cm.txt MAP2bp.txt \
#     2_analysis/ibd/no_outgr_direct_8.segments.tsv \
#     2_analysis/ibd/no_outgr_direct_8.ibd.tsv
# ===============================================================
# This script produces converts truffle segments based on the
# hamlet linkage maps
# ---------------------------------------------------------------
# ===============================================================
# args <- c("ressources/recombination/",
#           "MAP1cm.txt", "MAP1bp.txt",
#           "MAP2cm.txt", "MAP2bp.txt",
#           "2_analysis/ibd/no_outgr_direct_8.segments.tsv",
#           "2_analysis/ibd/no_outgr_direct_8.ibd.tsv")
# script_name <- "R/recombination_scaling.R"
args <- commandArgs(trailingOnly = FALSE)
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(hypogen)
library(patchwork)
library(plyranges)

cat('\n')
script_name <- args[5] %>%
  str_remove(., '--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ', getwd(), '/', .)

args <- process_input(script_name, args)
# config -----------------------
rec_path <- as.character(args[1])
hypo_map1_cm <- as.character(args[2])
hypo_map1_bp <- as.character(args[3])
hypo_map2_cm <- as.character(args[4])
hypo_map2_bp <- as.character(args[5])
segment_file <- as.character(args[6])
summary_file <- as.character(args[7])

truffle_idx <- segment_file %>% 
  str_extract("[0-9]*.segments.tsv") %>% 
  str_remove(".segments.tsv")
# custom functions ----------------
read_maps <- function(cm_file, bp_file){
  read_tsv(cm_file) %>% 
    group_by(LG) %>% 
    mutate(LGnm = as.roman(LG) %>% as.numeric(),
           CHROM = str_c("LG", str_pad(LGnm, width = 2, pad = 0)),
           cM = ifelse(LG == "VIII", max(cM)-cM,cM)) %>% 
    ungroup() %>% 
    dplyr::select(Loci, cM, CHROM)  %>% 
    full_join(read_tsv(bp_file) %>% 
                filter(!(duplicated(Loci) | duplicated(Loci, fromlast = TRUE))),
              by = c(Loci = "Loci", CHROM = "LG")) %>% 
    left_join(hypo_chrom_start) %>% 
    mutate(GPOS = GSTART + bp) %>% 
    filter(!is.na(GPOS),
           !is.na(cM)) %>% 
    arrange(GPOS)
}

make_lg_seg <- function(lg = "LG08", n = 31, gmap = gmap1){
  data_pos <- tibble(CHROM = rep(lg, n),
                     start = seq(from = hypo_karyotype$GSTART[hypo_karyotype$CHROM == lg],
                                 to = hypo_karyotype$GEND[hypo_karyotype$CHROM == lg],
                                 length = n) %>%
                       floor(),
                     GSTART = hypo_karyotype$GSTART[hypo_karyotype$CHROM == lg]) %>% 
    mutate(GPOS = start,
           start = start - GSTART,
           end = start) %>% 
    as_iranges()
  
  map_pos <- gmap %>%
    filter(CHROM == lg) %>% 
    attach_end(LG = lg) %>% 
    mutate(start = lag(bp, default = 0),
           start_cM = lag(cM, default = 0)) %>% 
    dplyr::select(CHROM, start, end = bp, start_cM, end_cM = cM) %>%
    mutate(start_bp = start, end_bp = end) %>% 
    as_iranges()
  
  list(data = data_pos, map = map_pos)
}

attach_end <- function(data, LG = "LG01"){
  data %>% 
    bind_rows(., 
              data %>%
                filter(row_number() == last(row_number())) %>% 
                mutate(GPOS = hypo_karyotype$GEND[hypo_karyotype$CHROM == CHROM],
                       bp = hypo_karyotype$LENGTH[hypo_karyotype$CHROM == CHROM],
                       Loci = as.numeric(
                         str_c("-99",
                               str_remove(string = CHROM, "LG"))
                       )
                ))
}

bin_rescaler <- function(bp, start_cM, end_cM, start_bp, end_bp,...){
  scales::rescale(x = bp,
                  to = c(start_cM, end_cM),
                  from = c(start_bp, end_bp))
}

interpol_data <- function(lg, ...){
  data_pair <- make_lg_seg(lg = lg, ...)
  
  plyranges::join_overlap_inner(data_pair$data,
                                data_pair$map) %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    dplyr::select(CHROM = CHROM.x, bp = start, GSTART, GPOS, start_cM:end_bp) %>% 
    mutate(interpol_cM = pmap_dbl(cur_data(), bin_rescaler)) 
}

na_to_zero <- function(x){x[is.na(x)] <- 0}

convert_bp_to_cm <- function(data, lg = "LG08", gmap = gmap1){
  gmap_in <- deparse(substitute(gmap))
  
  data_pos <- data %>% 
    filter( CHROM == lg ) %>% 
    dplyr::select(PAIR, TYPE, CHROM, START, END, NMARKERS) %>% 
    mutate(seg_id = str_c(PAIR,"_",CHROM,"_",START)) %>% 
    pivot_longer(cols = START:END, names_to = "PART", values_to = "start") %>%
    mutate(end = start) %>% 
    as_iranges()
  
  map_pos <- gmap %>%
    filter(CHROM == lg) %>%
    attach_end(LG = lg) %>%
    mutate(start = lag(bp, default = 0) + 1, # avoid overlapping segments - causes duplications in joining
           start_cM = lag(cM, default = 0)) %>%
    dplyr::select(start, end = bp, start_cM, end_cM = cM) %>%
    mutate(start_bp = start, end_bp = end) %>%
    as_iranges()
  
  map_nr <- str_replace(gmap_in, pattern = "gmap", replacement = "_m")
  
  plyranges::join_overlap_inner(data_pos,
                                map_pos) %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::select(CHROM = CHROM, bp = start, PAIR, TYPE, PART, start_cM:end_bp, seg_id, NMARKERS) %>%
    mutate(interpol_cM = pmap_dbl(cur_data(), bin_rescaler)) %>%
    pivot_wider(id_cols = c(CHROM,PAIR,TYPE,seg_id, PART,NMARKERS),
                values_from = c(bp,interpol_cM),
                names_from = PART) %>%
    dplyr::select(-seg_id) %>%
    mutate(length_bp = bp_END - bp_START,
           length_cM = interpol_cM_END - interpol_cM_START) %>% 
    set_names(value = c("CHROM", "PAIR", "TYPE", "NMARKERS", "bp_START", "bp_END",
                        str_c(c("interpol_cM_START", "interpol_cM_END"), map_nr),
                        "length_bp", str_c("length_cM", map_nr)))
}

# actual script -------------------
gmap1 <- read_maps(cm_file = str_c(rec_path, hypo_map1_cm),
                   bp_file = str_c(rec_path, hypo_map1_bp))
gmap2 <- read_maps(cm_file = str_c(rec_path, hypo_map2_cm),
                   bp_file = str_c(rec_path, hypo_map2_bp))

# p1 <- gmap1 %>% 
#   ggplot(aes(x = GPOS, y = cM, group = CHROM)) +
#   geom_hypo_LG() +
#   geom_line() +
#   scale_fill_hypo_LG_bg()+
#   scale_x_hypo_LG() +
#   theme_hypo()
# 
# p2 <- gmap2 %>% 
#   ggplot(aes(x = GPOS, y = cM, group = CHROM)) +
#   geom_hypo_LG() +
#   geom_line() +
#   scale_fill_hypo_LG_bg()+
#   scale_x_hypo_LG() +
#   theme_hypo()
# 
# p1 / p2

lgs <- 1:24 %>%
  str_pad(width = 2, pad = 0) %>%
  str_c("LG",.)

segments_individual_interpol_map1 <- lgs %>% map_dfr(interpol_data, n = 51)
segments_individual_interpol_map2 <- lgs %>% map_dfr(interpol_data, n = 51, gmap = gmap2)

# p1 <-  gmap1 %>% 
#   ggplot(aes(x = GPOS, y = cM)) +
#   geom_line(color = rgb(0,0,0,.2),aes(group = CHROM))+
#   geom_point(size = .4) +
#   geom_point(data = segments_individual_interpol_map1, 
#              aes(y = interpol_cM), 
#              color = "red",
#              shape = 1,
#              size = .8) +
#   scale_x_hypo_LG() +
#   theme_minimal()
# 
# p2 <-  gmap2 %>% 
#   ggplot(aes(x = GPOS, y = cM)) +
#   geom_line(color = rgb(0,0,0,.2),aes(group = CHROM))+
#   geom_point(size = .4) +
#   geom_point(data = segments_individual_interpol_map2, 
#              aes(y = interpol_cM), 
#                  color = "red",
#              shape = 1,
#              size = .8) +
#   scale_x_hypo_LG() +
#   theme_minimal()
# 
# p1 / p2

segments_individual <- vroom::vroom(segment_file) %>%
  mutate(LENGTH = LENGTH * 10^6,
         START = POS * 10^6,
         END = START + LENGTH,
         PAIR = str_c(ID1, "-", ID2))

segments_summary <- vroom::vroom(summary_file) %>%
  mutate(PAIR = str_c(ID1, "-", ID2))

# check if summary of individual segments reproduces the original truffle summary
# segments_individual %>%
#   group_by(PAIR, TYPE) %>%
#   summarise(seq_length = sum(LENGTH),
#             n_mark = sum(NMARKERS)) %>%
#   ungroup() %>%
#   pivot_wider(id_cols = PAIR, names_from = TYPE, values_from = seq_length:n_mark) %>% 
#   left_join(segments_summary) %>% 
#   mutate(n_mark_IBD2 = na_to_zero(n_mark_IBD2),
#          IBD0_manual = (NMARK - (n_mark_IBD1 + n_mark_IBD2)) / NMARK,
#          IBD1_manual = n_mark_IBD1 / NMARK,
#          IBD2_manual = n_mark_IBD2 / NMARK,
#          icheck_0 = IBD0_manual - IBD0,
#          icheck_1 = IBD1_manual - IBD1,
#          icheck_2 = IBD2 - IBD2)

bounds_gmap1 <- gmap1 %>% 
  group_by(CHROM) %>% 
  filter(cM == max(cM)) %>% 
  filter(bp == max(bp)) %>% 
  ungroup() %>% 
  dplyr::select(CHROM, cM) %>%
  mutate(GSTART_cM = cumsum(lag(cM, default = 0)),
         GEND_cM = cM + GSTART_cM,
         GMID_cM = (GSTART_cM + GEND_cM) / 2,
         grp = c("even", "odd")[ 1+row_number() %% 2 ])

bounds_gmap2 <- gmap2 %>% 
  group_by(CHROM) %>% 
  filter(cM == max(cM)) %>% 
  filter(bp == max(bp)) %>% 
  ungroup() %>% 
  dplyr::select(CHROM, cM) %>%
  mutate(GSTART_cM = cumsum(lag(cM, default = 0)),
         GEND_cM = cM + GSTART_cM,
         GMID_cM = (GSTART_cM + GEND_cM) / 2,
         grp = c("even", "odd")[ 1+row_number() %% 2 ])

hypo_all_starts <- hypo_karyotype %>% 
  dplyr::select(CHROM, GSTART, GEND) %>% 
  left_join(bounds_gmap1 %>% 
              dplyr::select(CHROM, GSTART_cM_m1 = GSTART_cM, GEND_cM_m1 = GEND_cM)) %>% 
  left_join(bounds_gmap2 %>% 
              dplyr::select(CHROM, GSTART_cM_m2 = GSTART_cM, GEND_cM_m2 = GEND_cM))

converted_segments <- lgs %>% 
  map_dfr(convert_bp_to_cm, data = segments_individual) %>% 
  left_join( lgs %>% 
               map_dfr(convert_bp_to_cm, data = segments_individual, gmap = gmap2) ) %>%
  left_join(hypo_all_starts) %>% 
  mutate(G_SEG_START = GSTART + bp_START,
         G_SEG_END = GSTART + bp_END,
         G_SEG_START_cM_m1 = GSTART_cM_m1 + interpol_cM_START_m1,
         G_SEG_END_cM_m1 = GSTART_cM_m1 + interpol_cM_END_m1,
         G_SEG_START_cM_m2 = GSTART_cM_m2 + interpol_cM_START_m2,
         G_SEG_END_cM_m2 = GSTART_cM_m2 + interpol_cM_END_m2)

# converted_segments %>% 
#   pivot_longer(cols = c(length_cM_m1, length_cM_m2),names_to = "map", values_to = "length_cM") %>% 
#   mutate(map = map %>% str_remove("length_cM_")) %>% 
#   ggplot() +
#   geom_hex(aes(x = length_bp * 10^-6, y = length_cM),
#            color = rgb(0,0,0,.2)) +
#   scale_fill_viridis_c(option = "C") +
#   facet_grid(map ~ TYPE) +
#   theme_minimal() +
#   theme(panel.background = element_rect(fill ="transparent", color = rgb(0,0,0,.2)))

# ggsave("~/Desktop/ibd_converted_lengths.pdf", width = 5, height = 4, device = cairo_pdf)
  
# converted_segments %>%
#   ggplot() +
#   geom_segment(aes(x = bp_START* 10^-6, xend = bp_END* 10^-6,
#                    y = interpol_cM_START_m1, yend = interpol_cM_END_m1, group = PAIR, color = "map1"),
#                arrow = arrow(type = "closed", length = unit(4,"pt")), alpha = .3) +
#   geom_segment(aes(x = bp_START* 10^-6, xend = bp_END* 10^-6,
#                    y = interpol_cM_START_m2, yend = interpol_cM_END_m2, group = PAIR, color = "map2"),
#                arrow = arrow(type = "closed", length = unit(4,"pt")), alpha = .3) +
#   labs(x = "position (Mb)", y = "position (cM)") +
#   facet_wrap(CHROM ~ .) +
#   theme_minimal()

# ggsave("~/Desktop/ibd_converted_positions.pdf", width = 12,height = 8, device = cairo_pdf)

# p_dens <- converted_segments %>% 
#   pivot_longer(cols = c(length_cM_m1, length_cM_m2),names_to = "map", values_to = "length_cM") %>% 
#   mutate(map = map %>% str_remove("length_cM_"))%>%
#   ggplot(aes(x = abs(length_cM))) +
#   geom_density(fill = rgb(0,0,0, .2)) +
#   geom_vline(xintercept = .1, color = "red", linetype = 3) +
#   geom_vline(xintercept = .165, linetype = 3) +
#   facet_grid(map ~.)
# 
# p_dens + 
#   p_dens +
#   coord_cartesian(xlim = c(0,3)) +
#   annotation_custom(grob = grid::textGrob(label = "McGee thresshold (0.1 cM)", gp = grid::gpar(col = "red")),ymin = .3) +
#   plot_layout(widths = c(.5,1)) &
#   theme_minimal()

# ggsave("~/Desktop/ibd_converted_thresholds.pdf", width = 8,height = 4, device = cairo_pdf)

# p_map1 <- bounds_gmap1 %>% 
#   ggplot() +
#   geom_rect(aes(xmin = GSTART_cM, xmax = GEND_cM, ymin = -Inf, ymax = Inf, fill = grp)) +
#   coord_cartesian(ylim = 0:1) +
#   scale_x_continuous(breaks = bounds_gmap1$GMID_cM,
#                      labels = bounds_gmap1$CHROM %>%
#                        str_remove("LG"),
#                      position = "top",
#                      limits = c(0, max(bounds_gmap1$GEND_cM)),expand = c(0, 0)) +
#   scale_fill_manual(values = c(odd = "transparent", even = rgb(.5,.5,.5,.3))) +
#   theme_minimal()
# 
# p_map2 <- bounds_gmap2 %>% 
#   ggplot() +
#   geom_rect(aes(xmin = GSTART_cM, xmax = GEND_cM, ymin = -Inf, ymax = Inf, fill = grp)) +
#   coord_cartesian(ylim = 0:1) +
#   scale_x_continuous(breaks = bounds_gmap2$GMID_cM,
#                      labels = bounds_gmap2$CHROM %>%
#                        str_remove("LG"),
#                      position = "top",
#                      limits = c(0, max(bounds_gmap2$GEND_cM)),expand = c(0, 0)) +
#   scale_fill_manual(values = c(odd = "transparent", even = rgb(.5,.5,.5,.3))) +
#   theme_minimal()
# 
# p_map1 / p_map2
# 
# (p_map1 +
#     geom_segment(data = converted_segments %>%
#                    filter(PAIR == "16_21-30nigpan-18151nigbel", length_cM_m1 > .1 | length_cM_m2 > .1),
#                  aes(x = G_SEG_START_cM_m1, xend = G_SEG_END_cM_m1, y = .5, yend = .5, color = "map1"), size = 5)) /
# (p_map2 +
#   geom_segment(data = converted_segments %>%
#                  filter(PAIR == "16_21-30nigpan-18151nigbel", length_cM_m1 > .1 | length_cM_m2 > .1),
#                aes(x = G_SEG_START_cM_m2, xend = G_SEG_END_cM_m2, y = .5, yend = .5, color = "map2"), size = 5)) &
#   scale_color_manual(values = c(map1 = "red", map2 = "blue"))

hypo_cM_length_map1 <- max(hypo_all_starts$GEND_cM_m1)
hypo_cM_length_map2 <- max(hypo_all_starts$GEND_cM_m2)
hypo_bp_length <- hypo_karyotype$GEND[hypo_karyotype$CHROM == "LG24"]

control <- converted_segments %>%
  ungroup() %>%
  group_by(PAIR, TYPE) %>%
  summarise(seq_length = sum(length_bp),
            n_mark = sum(NMARKERS),
            cm_length_m1 = sum(length_cM_m1),
            cm_length_m2 = sum(length_cM_m2)) %>%
  ungroup() %>%
  pivot_wider(id_cols = PAIR, names_from = TYPE, values_from = seq_length:cm_length_m2) %>% 
  left_join(segments_summary, .) %>% 
  mutate(n_mark_IBD2 = na_to_zero(n_mark_IBD2),
         IBD0_manual = (NMARK - (n_mark_IBD1 + n_mark_IBD2)) / NMARK,
         IBD1_manual = n_mark_IBD1 / NMARK,
         IBD2_manual = n_mark_IBD2 / NMARK,
         icheck_0 = IBD0_manual - IBD0,
         icheck_1 = IBD1_manual - IBD1,
         icheck_2 = IBD2 - IBD2,
         # compile ibd by sequence map
         ibd0_bp = (hypo_bp_length - (seq_length_IBD1 + seq_length_IBD2)) / hypo_bp_length,
         ibd1_bp = seq_length_IBD1 / hypo_bp_length,
         ibd2_bp = seq_length_IBD2 / hypo_bp_length,
         # compile ibd by genetic map 1
         ibd0_cM_m1 = (hypo_cM_length_map1 - (cm_length_m1_IBD1 + cm_length_m1_IBD2)) / hypo_cM_length_map1,
         ibd1_cM_m1 = cm_length_m1_IBD1 / hypo_cM_length_map1,
         ibd2_cM_m1 = cm_length_m1_IBD2 / hypo_cM_length_map1,
         # compile ibd by genetic map 2
         ibd0_cM_m2 = (hypo_cM_length_map2 - (cm_length_m2_IBD1 + cm_length_m2_IBD2)) / hypo_cM_length_map2,
         ibd1_cM_m2 = cm_length_m2_IBD1 / hypo_cM_length_map2,
         ibd2_cM_m2 = cm_length_m2_IBD2 / hypo_cM_length_map2)

cM_treshold <- 0.2
summary_filterd <- converted_segments %>%
  filter(length_cM_m1 > cM_treshold & length_cM_m2 > cM_treshold) %>% 
  ungroup() %>%
  group_by(PAIR, TYPE) %>%
  summarise(seq_length = sum(length_bp),
            n_mark = sum(NMARKERS),
            cm_length_m1 = sum(length_cM_m1),
            cm_length_m2 = sum(length_cM_m2)) %>%
  ungroup() %>%
  pivot_wider(id_cols = PAIR, names_from = TYPE, values_from = seq_length:cm_length_m2) %>% 
  left_join(segments_summary, .) %>% 
  mutate(n_mark_IBD2 = na_to_zero(n_mark_IBD2),
         IBD0_manual = (NMARK - (n_mark_IBD1 + n_mark_IBD2)) / NMARK,
         IBD1_manual = n_mark_IBD1 / NMARK,
         IBD2_manual = n_mark_IBD2 / NMARK,
         icheck_0 = IBD0_manual - IBD0,
         icheck_1 = IBD1_manual - IBD1,
         icheck_2 = IBD2 - IBD2,
         # compile ibd by sequence map
         ibd0_bp = (hypo_bp_length - (seq_length_IBD1 + seq_length_IBD2)) / hypo_bp_length,
         ibd1_bp = seq_length_IBD1 / hypo_bp_length,
         ibd2_bp = seq_length_IBD2 / hypo_bp_length,
         # compile ibd by genetic map 1
         ibd0_cM_m1 = (hypo_cM_length_map1 - (cm_length_m1_IBD1 + cm_length_m1_IBD2)) / hypo_cM_length_map1,
         ibd1_cM_m1 = cm_length_m1_IBD1 / hypo_cM_length_map1,
         ibd2_cM_m1 = cm_length_m1_IBD2 / hypo_cM_length_map1,
         # compile ibd by genetic map 2
         ibd0_cM_m2 = (hypo_cM_length_map2 - (cm_length_m2_IBD1 + cm_length_m2_IBD2)) / hypo_cM_length_map2,
         ibd1_cM_m2 = cm_length_m2_IBD1 / hypo_cM_length_map2,
         ibd2_cM_m2 = cm_length_m2_IBD2 / hypo_cM_length_map2)

control_subplot <- function(x,y,c, data){
    ggplot(data) +
    geom_point(aes_string(x = x, y = y, color = str_c("'", c, "'")),
               alpha = .4)
}

control_subplot("IBD0", "IBD0_manual", "IBD0") +
  control_subplot("IBD1", "IBD1_manual", "IBD1") 

replacer <-  function(str, ibd){str_replace(str, pattern = "([IBDibd]{3})X", replacement = str_c("\\1",ibd))}
control_plot <- function(x, y, cl, data){
  
  xs <- map_chr(0:2, .f = replacer, str = x)
  ys <- map_chr(0:2, .f = replacer, str = y)
  cls <- map_chr(0:2, .f = replacer, str = cl)
  
  clrs <- RColorBrewer::brewer.pal(3, "Set1") %>%
    set_names(value = cls)

  control_subplot(xs[1], ys[1], cls[1], data = data) +
    control_subplot(xs[2], ys[2], cls[2], data = data) +
    control_subplot(xs[3], ys[3], cls[3], data = data) +
    plot_layout(nrow = 1, 
                guides = "collect") &
    scale_color_manual(values = clrs, guide = "none") &
    theme_minimal()
}

to_check <- c("IBDX", "IBDX_manual", "ibdX_bp", "ibdX_cM_m1", "ibdX_cM_m2")


control_plot(to_check[1], to_check[2], to_check[2], data = control)/
control_plot(to_check[1], to_check[3], to_check[3], data = control)/
control_plot(to_check[1], to_check[4], to_check[4], data = control)/
control_plot(to_check[1], to_check[5], to_check[5], data = control) 

ggsave("~/Desktop/ibd_comparisons.pdf", width = 8, height = 8, device = cairo_pdf)
ggsave("~/Desktop/ibd_comparisons.png", width = 8, height = 8, type = "cairo")

control_plot(to_check[1], to_check[2], to_check[2], data = summary_filterd)/
  control_plot(to_check[1], to_check[3], to_check[3], data = summary_filterd)/
  control_plot(to_check[1], to_check[4], to_check[4], data = summary_filterd)/
  control_plot(to_check[1], to_check[5], to_check[5], data = summary_filterd) +
  plot_annotation(title = glue::glue("IBD segments filtered by length > {cM_treshold}"))

ggsave("~/Desktop/ibd_comparisons_filtered.png", width = 8, height = 8, type = "cairo")
