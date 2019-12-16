#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SFx9.R 2_analysis/het/50kb/
# ===============================================================
# by Hench, McMillan an Puebla
#   ---------------------------------------------------------------
# ===============================================================
# args <- c('2_analysis/het/50kb/')
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(hypoimg)
library(GenomicOriginsScripts)


cat('\n')
script_name <- args[5] %>%
  str_remove(.,'--file=')

plot_comment <- script_name %>%
  str_c('mother-script = ',getwd(),'/',.)

args <- process_input(script_name, args)

# config -----------------------
het_dir <- as.character(args[1])

files <- dir(het_dir)

collapse_bins <- function(tib){
  tib %>%
    mutate(win_id = str_c(CHROM,BIN_START, BIN_END, sep = '_')) %>%
    select(-CHROM,-BIN_START, -BIN_END, -N_SNPs)
}
#read_tsv(str_c(het_dir,'20644abehon.50k.5k.het.tsv.gz'))
data <- str_c(het_dir, files) %>%
  map(vroom, delim = '\t') %>%
  map(collapse_bins) %>%
  reduce(left_join) %>% 
  pivot_longer(names_to = 'ind', values_to = 'het', cols = -win_id) %>%
  separate(win_id, 
           into = c('CHROM', 'win_start', 'win_end'),
           sep = '_',convert = TRUE) %>%
  left_join(hypo_karyotype) %>%
  mutate(gpos = GSTART + (win_start + win_end)/2,
         spec = str_sub(ind,-6,-4),
         loc = str_sub(ind,-3,-1))

grp_summary <- data %>%
  group_by(CHROM, gpos, spec, loc) %>%
  summarise(avg_het = mean(het, na.rm = TRUE),
            min_het = min(het, na.rm = TRUE),
            max_het = max(het, na.rm = TRUE),
            sd_het = sd(het, na.rm = TRUE)) %>%
  ungroup()

clr_alt <- clr
clr_alt['uni'] <- rgb(.7,.7,.7)

p_het_wg <- data %>% 
  ggplot() +
  geom_hypo_LG() +
  geom_line(mapping = aes(x = gpos, y = het, group = ind,
                          color = spec), alpha = .4)+
  scale_x_hypo_LG()+
  scale_y_continuous(name = 'Heterozygosity', 
                     breaks = c(0,.01,.02))+
  scale_fill_hypo_LG_bg()+
  scale_color_manual(values = clr_alt)+
  facet_wrap(facets = vars(spec, loc),
             strip.position = 'right', ncol = 1)+
  #coord_cartesian(ylim = c(0,.5))+
  theme_hypo()+
  theme(legend.position = 'none',
        strip.background = element_blank())

hypo_save(plot = p_het_wg,
          filename = 'figures/SX9a.png',
          width = 12,
          height = 15,
          dpi = 200,
          comment = plot_comment)

p_het_lg06 <- data %>% 
  filter(CHROM == 'LG06') %>%
  ggplot() +
  geom_line(mapping = aes(x = gpos, y = het, group = ind,
                          color = spec), alpha = .4)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(name = 'Heterozygosity', 
                     breaks = c(0,.01,.02))+
  scale_fill_hypo_LG_bg()+
  scale_color_manual(values = clr_alt)+
  facet_wrap(facets = vars(spec, loc),
             strip.position = 'right', ncol = 1)+
  #coord_cartesian(ylim = c(0,.5))+
  theme_hypo()+
  theme(legend.position = 'none',
        strip.background = element_blank())

hypo_save(plot = p_het_lg06,
          filename = 'figures/SX9b.png',
          width = 10,
          height = 10,
          dpi = 200,
          comment = plot_comment)
  
# 
# grp_summary %>% 
#   filter(!(spec %in% c('tor', 'tab', 'flo'))) %>%
#   group_by(spec) %>%
#   summarise(genome_wide_het = mean(avg_het)) %>%
#   arrange(genome_wide_het)
