#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF3.R \
#   2_analysis/summaries/fst_globals.txt \
#   2_analysis/fst/50k/ \
#   2_analysis/fasteprr/step4/fasteprr.all.rho.txt.gz
# ===============================================================
# This script
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "2_analysis/admixture/", "metadata/phenotypes.sc")
args <- commandArgs(trailingOnly=FALSE)
# setup -----------------------
library(paletteer)
#library(ggstance)
library(patchwork)
library(GenomicOriginsScripts)
library(hypoimg)
library(ggtext)

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
admx_path <- as.character(args[1])
pheno_file <- as.character(args[2])


gids <- dir(admx_path,
            pattern = "pop.*15.txt") %>% 
  str_remove("pop.") %>% 
  str_remove(".15.txt") 

pheno_data <- read_sc(pheno_file) %>%
  select(id, Bars, Peduncle, Snout ) %>%
  filter(!is.na(Bars)) 


data <- gids %>%
  map_dfr(data_amdx,admx_path = admx_path,
          k = 2) 

pheno_facet <- tibble( trait = c("Snout","Bars",  "Peduncle"),
                       gid = c("LG04_1", "LG12_3", "LG12_4")) %>%
  mutate(facet_label = str_c(gid, " / ", trait))


gid_labels <-  c(LG04_1 = "LG04 (A)",
                 LG12_3 = "LG12 (B)", 
                 LG12_4 = "LG12 (C)")

gid_traits <-  c(LG04_1 = "Snout",
                 LG12_3 = "Bars", 
                 LG12_4 = "Peduncle")

pheno_plot_data <- data %>%
  filter(!duplicated(id)) %>%
  select(id:id_order) %>%
  left_join(pheno_data,by = c( id_nr = "id")) %>%
  arrange(spec, Bars, Peduncle, Snout, id) %>%
  mutate(ord_nr = row_number()) %>% 
  pivot_longer(names_to = "trait",
               values_to = "phenotype",
               cols = Bars:Snout) %>%
  left_join(pheno_facet)

sample_order <- pheno_plot_data %>%
  filter(!duplicated(id)) %>%
  select(id, ord_nr)

p_ad <- c("LG04_1", "LG12_3", "LG12_4") %>% map(adm_plot)

p_phno <- pheno_plot_data %>%
  ggplot(aes(x = factor(ord_nr)))+
  geom_point(aes(y = trait, fill = factor(phenotype)),shape = 21)+
  # scale_y_continuous(expand = c(0,0))+#,limits = c(-.1,1))+
  scale_x_discrete(breaks = sample_order$ord_nr,
                   labels = sample_order$id) +
  scale_fill_manual("Phenotype<br><img src='ressources/img/all_traits_c.png' width='110' />",
                    values = c(`0` = "white", `1` = "black"),
                    na.value = "gray",
                    labels = c("absent", "present", "not scored"))+
  guides(fill = guide_legend(ncol = 1))+
  theme_minimal()+
  theme(#axis.text.y = element_text(size = 4),
    plot.title = element_text(size = 9),
    legend.title = element_markdown(hjust = .5),
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank()#element_text(angle = 90, size = 6)
  )

tib_drawing <- pheno_plot_data %>%
  group_by(spec) %>%
  summarise(pos = (min(ord_nr)+max(ord_nr))*.5) %>%
  ungroup() 

p_spec <- pheno_plot_data %>%
  group_by(spec) %>%
  summarise(start = min(ord_nr)-1,
            end = max(ord_nr)) %>%
  ggplot(aes(xmin = start, xmax = end,
             ymin = -Inf, 
             ymax = Inf))+
  geom_rect(aes(fill = spec), color = "black")+
  (tib_drawing %>% pmap(add_spec_drawing))+
  scale_y_continuous(breaks = .5,labels = c( "Species"), limits = c(0,1))+#,limits = c(-.1,1))+
  scale_x_discrete(breaks = sample_order$ord_nr,
                   labels = sample_order$id,
                   expand = c(0,0)) +
  scale_fill_manual("Species", values = clr, labels = sp_labs)+
  theme_minimal()+
  theme(#axis.text.y = element_text(size = 4),
    plot.title = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank()#element_text(angle = 90, size = 6)
  )

p_loc <- pheno_plot_data %>%
  ggplot(aes(x = factor(ord_nr)))+
  geom_raster(aes(y = 0, fill = loc))+
  scale_y_continuous(breaks = c(0),labels = c("Location"))+#,limits = c(-.1,1))+
  scale_x_discrete(breaks = sample_order$ord_nr,
                   labels = sample_order$id) +
  scale_fill_manual("Location", values =  clr_loc, loc_names)+
  theme_minimal()+
  theme(#axis.text.y = element_text(size = 4),
    plot.title = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank()#element_text(angle = 90, size = 6)
  )


p_l <- (get_legend(p_phno) %>% ggdraw()) +
  (get_legend(p_spec) %>% ggdraw()) +
  (get_legend(p_loc) %>% ggdraw()) +
  plot_layout(nrow = 1)

p_done <- p_ad[[1]] + 
  p_ad[[2]] +
  p_ad[[3]]+ 
  p_spec + p_loc + p_l + 
  plot_layout(ncol = 1, heights = c(.4,.4,.4,.08,.02,.1)) & 
  theme(legend.position = "none")

scl <- .9
ggsave("figures/SF10.pdf",
       plot = p_done,
       width = 16*scl,
       height = 10*scl,
       device = cairo_pdf)

