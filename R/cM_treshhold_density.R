args <- c("ressources/recombination/",
          "MAP1cm.txt", "MAP1bp.txt",
          "MAP2cm.txt", "MAP2bp.txt",
          "10")
# setup -----------------------
renv::activate()
library(GenomicOriginsScripts)
library(hypogen)
library(patchwork)
library(plyranges)

rec_path <- as.character(args[1])
hypo_map1_cm <- as.character(args[2])
hypo_map1_bp <- as.character(args[3])
hypo_map2_cm <- as.character(args[4])
hypo_map2_bp <- as.character(args[5])
idx <- as.character(as.character(args[6]))

segment_file <- glue::glue("2_analysis/ibd/no_outgr_bed95_{idx}.segments.tsv")
summary_file <- glue::glue("2_analysis/ibd/no_outgr_bed95_{idx}.ibd.tsv")

iterations <- c(str_c("2/5*10^",6:3," BP"),"25/10 kb","10/5 kb","7-5/3 kb","15/7.5 kb")
itteration_names <- c(str_c("10-",6:3),"7","8","9","10")

tag <- iterations[which(itteration_names == idx)] %>% str_replace("/", "_") %>% str_replace_all("[ \\.]","-")

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

na_to_zero <- function(x){
  x_type <- typeof(x)
  if_else(is.na(x), as(0,Class = x_type), x) %>% 
    as.double() %>% as(Class = x_type)
  
}

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


lgs <- 1:24 %>%
  str_pad(width = 2, pad = 0) %>%
  str_c("LG",.)

segments_individual_interpol_map1 <- lgs %>% map_dfr(interpol_data, n = 51)
segments_individual_interpol_map2 <- lgs %>% map_dfr(interpol_data, n = 51, gmap = gmap2)

segments_individual <- vroom::vroom(segment_file) %>%
  mutate(LENGTH = LENGTH * 10^6,
         START = POS * 10^6,
         END = START + LENGTH,
         PAIR = str_c(ID1, "-", ID2))

segments_summary <- vroom::vroom(summary_file) %>%
  mutate(PAIR = str_c(ID1, "-", ID2))

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

p_dens <- converted_segments %>%
  pivot_longer(cols = c(length_cM_m1, length_cM_m2),
               names_to = "map", values_to = "length_cM") %>%
  mutate(map = map %>% str_remove("length_cM_")) %>%
  ggplot(aes(x = abs(length_cM))) +
  geom_density(fill = rgb(0,0,0, .2)) +
  # geom_vline(xintercept = .1, color = "red", linetype = 3) +
  geom_vline(xintercept = .2, linetype = 3, color = "red") +
  facet_grid(map ~.)

p_dens +
  p_dens +
  coord_cartesian(xlim = c(0,3)) +
  plot_layout(widths = c(.5,1)) &
  theme_minimal()

ggsave(glue::glue("~/Desktop/ibd_converted_thresholds_{tag}.pdf"), width = 8, height = 4, device = cairo_pdf)
