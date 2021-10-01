library(GenomicOriginsScripts)
library(hypogen)
library(patchwork)
library(plyranges)

read_maps <- function(cm_file, bp_file){
  read_tsv(cm_file) %>% 
    group_by(LG) %>% 
    mutate(LGnm = as.roman(LG) %>% as.numeric(),
           CHROM = str_c("LG", str_pad(LGnm, width = 2, pad = 0)),
           cM = ifelse(LG == "VIII", max(cM)-cM,cM)) %>% 
    ungroup() %>% 
    dplyr::select(Loci, cM, CHROM)  %>% 
    full_join(read_tsv(bp_file) %>% 
                filter(!(duplicated(Loci) | duplicated(Loci,fromlast = TRUE))),
              by = c(Loci = "Loci", CHROM = "LG")) %>% 
    left_join(hypo_chrom_start) %>% 
    mutate(GPOS = GSTART + bp) %>% 
    filter(!is.na(GPOS),
           !is.na(cM)) %>% 
    arrange(GPOS)
}

rec_path <- "ressources/recombination/"
gmap1 <- read_maps(cm_file = str_c(rec_path, "MAP1cm.txt"), bp_file = str_c(rec_path, "MAP1bp.txt"))
gmap2 <- read_maps(cm_file = str_c(rec_path, "MAP2cm.txt"), bp_file = str_c(rec_path, "MAP2bp.txt"))

p1 <- gmap1 %>% 
  ggplot(aes(x = GPOS, y = cM, group = CHROM)) +
  geom_hypo_LG() +
  geom_line() +
  scale_fill_hypo_LG_bg()+
  scale_x_hypo_LG() +
  theme_hypo()

p2 <- gmap2 %>% 
  ggplot(aes(x = GPOS, y = cM, group = CHROM)) +
  geom_hypo_LG() +
  geom_line() +
  scale_fill_hypo_LG_bg()+
  scale_x_hypo_LG() +
  theme_hypo()

p1 / p2


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

lgs <- 1:24 %>%
  str_pad(width = 2, pad = 0) %>%
  str_c("LG",.)

test_interpol_map1 <- lgs %>% map_dfr(interpol_data, n = 51)
test_interpol_map2 <- lgs %>% map_dfr(interpol_data, n = 51, gmap = gmap2)

p1 <-  gmap1 %>% 
  ggplot(aes(x = GPOS, y = cM)) +
  geom_line(color = rgb(0,0,0,.2),aes(group = CHROM))+
  geom_point(size = .4) +
  geom_point(data = test_interpol_map1, 
             aes(y = interpol_cM), 
             color = "red",
             shape = 1,
             size = .8) +
  scale_x_hypo_LG() +
  theme_minimal()

p2 <-  gmap2 %>% 
  ggplot(aes(x = GPOS, y = cM)) +
  geom_line(color = rgb(0,0,0,.2),aes(group = CHROM))+
  geom_point(size = .4) +
  geom_point(data = test_interpol_map2, 
             aes(y = interpol_cM), 
                 color = "red",
             shape = 1,
             size = .8) +
  scale_x_hypo_LG() +
  theme_minimal()

p1 / p2

# --------- #

na_to_zero <- function(x){x[is.na(x)] <- 0}



test <- vroom::vroom("2_analysis/ibd/no_outgr_direct_8.segments.tsv") %>%
  mutate(LENGTH = LENGTH * 10^6,
         START = POS * 10^6,
         END = START + LENGTH,
         PAIR = str_c(ID1, "-", ID2))

test_2 <- vroom::vroom("2_analysis/ibd/no_outgr_direct_8.ibd.tsv") %>%
  mutate(PAIR = str_c(ID1, "-", ID2))

convert_bp_to_cm <- function(data, lg = "LG08", gmap = gmap1){
  gmap_in <- deparse(substitute(gmap))
  data_pos <- data %>% 
    filter( CHROM == lg) %>% 
    dplyr::select(PAIR, TYPE, CHROM, START, END, NMARKERS) %>% 
    mutate(seg_id = str_c(PAIR,"_",CHROM,"_",START)) %>% 
    pivot_longer(cols = START:END, names_to = "PART", values_to = "start") %>%
    mutate(end = start) %>% 
    as_iranges()
  
  map_pos <- gmap %>%
    filter(CHROM == lg) %>% 
    attach_end(LG = lg) %>% 
    mutate(start = lag(bp, default = 0),
           start_cM = lag(cM, default = 0)) %>% 
    dplyr::select(start, end = bp, start_cM, end_cM = cM) %>%
    mutate(start_bp = start, end_bp = end) %>% 
    as_iranges()
  
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
           length_cM = interpol_cM_END - interpol_cM_START,
           map = gmap_in)
}

converted_segments <- lgs %>% 
  map_dfr(convert_bp_to_cm, data = test) %>% 
  bind_rows( lgs %>% 
               map_dfr(convert_bp_to_cm, data = test, gmap = gmap2))

converted_segments %>% 
  ggplot() +
  geom_hex(aes(x = length_bp * 10^-6, y = length_cM),
           color = rgb(0,0,0,.2)) +
  scale_fill_viridis_c(option = "C") +
  facet_grid(map ~ TYPE) +
  theme_minimal() +
  theme(panel.background = element_rect(fill ="transparent", color = rgb(0,0,0,.2)))

ggsave("~/Desktop/ibd_converted_lengths.pdf", width = 5, height = 4, device = cairo_pdf)
  
converted_segments %>%
  ggplot() +
  geom_segment(aes(x = bp_START* 10^-6, xend = bp_END* 10^-6,
                   y = interpol_cM_START, yend = interpol_cM_END, group = PAIR, color = map),
               arrow = arrow(type = "closed", length = unit(4,"pt")), alpha = .3) +
  labs(x = "position (Mb)", y = "position (cM)") +
  facet_wrap(CHROM ~ .) +
  theme_minimal()

ggsave("~/Desktop/ibd_converted_positions.pdf", width = 12,height = 8, device = cairo_pdf)

p_dens <- converted_segments %>%
  ggplot(aes(x = abs(length_cM))) +
  geom_density(fill = rgb(0,0,0, .2)) +
  geom_vline(xintercept = .1, color = "red", linetype = 3) +facet_grid(map ~.)

p_dens + 
  p_dens +
  coord_cartesian(xlim = c(0,3)) +
  annotation_custom(grob = grid::textGrob(label = "McGee thresshold (0.1 cM)", gp = grid::gpar(col = "red")),ymin = .3) +
  plot_layout(widths = c(.5,1)) &
  theme_minimal()
ggsave("~/Desktop/ibd_converted_thresholds.pdf", width = 8,height = 4, device = cairo_pdf)


test %>%
  group_by(PAIR, TYPE) %>%
  summarise(seq_length = sum(LENGTH),
            n_mark = sum(NMARKERS)) %>%
  ungroup() %>%
  pivot_wider(id_cols = PAIR, names_from = TYPE, values_from = seq_length:n_mark) %>% 
  left_join(test_2) %>% 
  mutate(n_mark_IBD2 = na_to_zero(n_mark_IBD2),
         IBD0_manual = (NMARK - (n_mark_IBD1 + n_mark_IBD2)) / NMARK,
         IBD1_manual = n_mark_IBD1 / NMARK,
         IBD2_manual = n_mark_IBD2 / NMARK,
         icheck_0 = IBD0_manual - IBD0,
         icheck_1 = IBD1_manual - IBD1,
         icheck_2 = IBD2 - IBD2)

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

p_map1 <- bounds_gmap1 %>% 
  ggplot() +
  geom_rect(aes(xmin = GSTART_cM, xmax = GEND_cM, ymin = -Inf, ymax = Inf, fill = grp)) +
  coord_cartesian(ylim = 0:1) +
  scale_x_continuous(breaks = bounds_gmap1$GMID_cM,
                     labels = bounds_gmap1$CHROM %>%
                       str_remove("LG"),
                     position = "top",
                     limits = c(0, max(bounds_gmap1$GEND_cM)),expand = c(0, 0)) +
  scale_fill_manual(values = c(odd = "transparent", even = rgb(.5,.5,.5,.3))) +
  theme_minimal()

p_map2 <- bounds_gmap2 %>% 
  ggplot() +
  geom_rect(aes(xmin = GSTART_cM, xmax = GEND_cM, ymin = -Inf, ymax = Inf, fill = grp)) +
  coord_cartesian(ylim = 0:1) +
  scale_x_continuous(breaks = bounds_gmap2$GMID_cM,
                     labels = bounds_gmap2$CHROM %>%
                       str_remove("LG"),
                     position = "top",
                     limits = c(0, max(bounds_gmap2$GEND_cM)),expand = c(0, 0)) +
  scale_fill_manual(values = c(odd = "transparent", even = rgb(.5,.5,.5,.3))) +
  theme_minimal()

p_map1 / p_map2
