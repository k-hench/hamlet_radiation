# -------------------------------------------------

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