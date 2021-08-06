library(tidyverse)
library(tidygraph)
library(ggraph)
library(prismatic)
library(patchwork)

library(fuzzyjoin)

n_zeros <- 4

iterations <- c(str_c("2/5*10^",6:3,"BP"),"25/10KB","10/5KB","7-5/3KB","15/7.5KB")
# itteration_names <-  c(as.character(6:3),"7","8","9","10")# 
itteration_names <-   c(str_c("10-",6:3),"7","8","9","10")

plot_network <- function(n_zeros, filt = 0){
  data <- read_tsv(glue::glue("output/no_outgr_{itteration_names[n_zeros]}.ibd.tsv")) %>% 
    mutate(ibd_total = (IBD2 + 0.5*IBD1) / (IBD0 + IBD1 + IBD2)) 
  
  data %>% 
    as_tbl_graph() %E>%
    filter(ibd_total > filt) %N>%
    mutate(spec = str_sub(name,-6,-4),
           loc = str_sub(name,-3,-1))  %>% 
    ggraph( layout = 'kk', weights = ibd_total) +
    geom_edge_link(aes(alpha = ibd_total), color = rgb(.1,.1,.1)) +
    geom_node_point(aes(color = spec,
                        shape = loc, fill = after_scale(clr_lighten(color,.2))), size = 3) +
    labs(title = glue::glue("seq length: {iterations[n_zeros]}, > {filt} IBD")) +
    scale_color_manual(values = GenomicOriginsScripts::clr) +
    scale_edge_alpha_continuous(limits = c(0,.3), guide = "none") +
    scale_shape_manual(values = 21:23) +
    coord_equal()
}

c(3,5,8,6,7,4) %>% map(plot_network) %>% wrap_plots(guides = "collect") 

ggsave("ibd_network_10-6-between_10-4_and_10-3_kk.png", width = 11, height = 8)  
ggsave("ibd_network_defaults.png", width = 5, height = 5)  

plot_hist <- function(n_zeros){
  data <- read_tsv(glue::glue("output/no_outgr_10-{n_zeros}.ibd.tsv")) %>% 
    mutate(ibd_total = (IBD2 + 0.5*IBD1) / (IBD0 + IBD1 + IBD2))  %>% 
  ggplot(aes(x = ibd_total)) +
  geom_histogram(bins = 40, fill = rgb(.7,.7,.7,.7), color = rgb(.7,.7,.7)) +
  labs(title = glue::glue("seq length: 2 & 5 * 10^{n_zeros} BP")) +
  theme_minimal()
}

6:3 %>% map(plot_hist) %>% wrap_plots(guides = "collect") 

ggsave("ibd_dist_10-6--10-3.pdf", width = 7,height = 4, device = cairo_pdf)

data_seg <- vroom::vroom("output/no_outgr_10.segments.tsv", delim = "\t") %>% 
  left_join(hypogen::hypo_chrom_start) %>% 
  mutate(start = POS * 10^6,
         end = start + (LENGTH * 10^6),
         gstart = GSTART + start,
         gend = start + (LENGTH * 10^6),
         ibd_hplo = str_remove(TYPE,"IBD") %>%
           as.integer())

library(IRanges)
library(plyranges)
library(hypogen)

hap_to_perc <- 100 / (166 * 165)

outlier_regions <- read_tsv("~/work/puebla_lab/chapter2_collect/hamlet_radiation/2_analysis/summaries/fst_outliers_998.tsv")

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
  geom_vline(data = outlier_regions, aes(xintercept = gpos), color = rgb(1,0,0,.4), size = .3) +
  geom_line(aes(x = GPOS, y = score * hap_to_perc, group = CHROM), color = rgb(.3,.3,.3), size = .3) +
  geom_ribbon(aes(x = GPOS, ymin = 0, ymax = score * hap_to_perc, group = CHROM), fill = rgb(0,0,0,.2)) +
  scale_hypobg_manual(values = c("transparent",rgb(.9,.9,.9,.9),"red","blue") %>%
                        set_names(nm = c("even", "odd", "a","b")), guide = "none")+
  scale_x_hypo_LG() +
  labs(y = "IBD score (percentage haplotypes in IBD)" 
         # "IBD score (coverage n haplotypes in IBD)"
       ) +
  coord_cartesian(ylim = c(-.001,4),#c(-2, 1000),
                  expand = 0) +
  theme_hypo()

ggsave("genome_wide_ibd_15_7-5kb.pdf", width = 13, height = 5, device = cairo_pdf)

