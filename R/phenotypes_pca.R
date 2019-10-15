#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla phenotypes_pca.R phenotypes.sc
# ===============================================================
# This script provides the pca based on the hamlet phenotypes
# ---------------------------------------------------------------
# ===============================================================
# args <- c('phenotypes.sc')
args = commandArgs(trailingOnly=FALSE)
args = args[7:8]
print(args)
# setup -----------------------
library(GenomicOriginsScripts)
library(logisticPCA)
library(ggrepel)

# config -----------------------
pheno_file <- as.character(args[1])

loading_scale <- 17
light_clr <- 'lightgray'
jit <- .4

# functions ------------
read_sc <-function(file){
  length_comment <- read_lines(file) %>% str_which(.,'^#') %>% max()
  length_format <- read_lines(file) %>% str_which(.,'^format*') %>% max()
  pos_frame <- read_lines(file) %>% str_which(.,'^frame*') %>% max()

  head_length <- max(length_comment, length_format, pos_frame)

  read_delim(file,delim = ' ',skip = head_length,
             col_names = c('type','cell','delim','value')) %>%
    filter(!(type == "goto"), !(is.na(value))) %>%
    select(-type,-delim) %>%
    mutate(cell = str_replace(cell,"^([A-Z])([0-9])","\\1_\\2")) %>%
    separate(cell, into = c('column','row'),sep = '_') %>%
    mutate(row = as.numeric(row)) %>%
    spread(key = column, value = value) %>%
    setNames(.,nm = .[1,] %>% as.character()) %>%
    filter(`0` != 0) %>%
    select(-`0`) %>%
    mutate_all(parse_guess)
}

inner <- function(data){data %>% select(Bars:Tail_transparent)}
transform_PCs <- function(x){scales::rescale(x = x, from = range(x),to = c(0,1))}
# the actual script -----------------

data <- read_sc(pheno_file) %>%
  select(-Tail_rim)

di <- data %>% inner()

# PCA (2d) ---------------
logsvd_model <- logisticSVD(di, k = 2)
logpca_cv <- cv.lpca(di, ks = 2, ms = 1:20)
logpca_model <- logisticPCA(di, k = 2, m = which.min(logpca_cv))

pca_scores <- logpca_model$PCs %>%
  as_tibble() %>%
  setNames(.,nm  = str_c('PC',1:2)) %>%
  bind_cols((data %>% select(id:geo)))

pca_lodings <- logpca_model$U %>%
  as_tibble() %>%
  setNames(.,nm  = str_c('PC',1:2)) %>%
  bind_cols(tibble(trait = names(di),
                   dir = c(1,1,1,1,-1,1,1,-1)))

# get convex hull of species on PCA
pca_scores_group <- pca_scores %>%
  select(spec,PC1,PC2) %>%
  group_by(spec) %>%
  slice(chull(PC1, PC2))
# get midpoint of species (for labels)
pca_scores_labs <- pca_scores_group %>%
  summarise(PC1 = mean(range(PC1)),
            PC2 = mean(range(PC2)))

plt_pc2 <- ggplot(pca_scores,aes(x=PC1,y=PC2))+
  geom_polygon(data = pca_scores_group,
               aes(fill=spec),alpha= .2,col=light_clr)+
  geom_jitter(aes(fill=spec),shape=21,size=2,width=jit,height=jit)+
  geom_segment(data = pca_lodings,inherit.aes = FALSE,
               aes(x=0,y=0,xend = PC1*loading_scale*dir,yend=PC2*loading_scale*dir),
               arrow = arrow(length = unit(3,'pt'),type = 'closed'),
               color=light_clr) +
  geom_text_repel(data = pca_lodings,inherit.aes = FALSE,segment.colour = NA,
                  aes(x = PC1*loading_scale*dir, y=PC2*loading_scale*dir,label=trait))+

  geom_text(data = pca_scores_labs,aes(label=spec),fontface='italic',color='red')+
  scale_fill_manual('Species',values = clr)+
  theme_minimal()+
  theme(panel.grid = element_blank(),axis.line = element_line(color = light_clr))

ggsave(plt_pc2, filename = 'phenotype_pca.pdf', width = 7, height = 6.5, device = cairo_pdf)

# PCA (1d) ---------------
logpca_model_D1 <- logisticPCA(di, k = 1, m = which.min(logpca_cv))
pca_scores_D1 <- logpca_model_D1$PCs %>%
  as_tibble() %>%
  setNames(.,nm  = str_c('PC',1)) %>%
  bind_cols((data %>% select(id:geo)))

labs <- pca_scores_D1 %>%
  group_by(spec) %>%
  summarise(mPC1 = mean(range(PC1)))

D1_data  <- pca_scores_D1 %>%
  left_join(labs)

p <- ggplot(D1_data,aes(x=PC1,y=1))+
  theme(axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

p1 <- p +
  geom_segment(aes(xend = mPC1,yend = 1.7,color = spec),alpha=.3)+
  geom_point(aes(fill=spec),shape=21,size=2)+
  geom_text_repel(inherit.aes = FALSE,data=labs,aes(y=2,x=mPC1,label=spec))+
  theme(legend.position = 'none',
        axis.title.x = element_blank())

p2 <- p + geom_jitter(aes(fill=spec),shape=21,size=2,width=0,height=jit)+
  theme(legend.position = 'bottom')

plt_pc1 <- cowplot::plot_grid(p1,p2,ncol = 1,align = 'v', rel_heights = c(.2,1))
ggsave(plt_pc1, filename = 'phenotype_pca_one_dimensional.pdf', width = 7, height = 8)

# export_table --------------------
# secies as phenotype
is_outgroup <- function(x){x$tab[row_number(x)] == 1 | x$tor[row_number(x)] == 1}
is_outgroup <- function(data,x,column){ifelse(!(data$tab[x] == 1 | data$tor[x]),column,NA)}
spec_pheno <- data %>%
  select(id,spec) %>%
  mutate(value = 1) %>%
  spread(value = value, key = spec, fill = 0) %>%
  mutate_at(.vars = vars(abe, gum, ind, may, nig, pue, ran, uni),
            .funs = function(x) ifelse(.$tor == 1 | .$tor == 1 , NA, x)) %>%
  select(-flo,-tab,-tor) # remove outgroups

data_export <- data %>%
  left_join( pca_scores %>% select(PC1:id) ) %>%
  left_join( pca_scores_D1  %>%
               select(PC1:id) %>%
               setNames(.,nm = c('PC_d1','id'))) %>%
  left_join(spec_pheno) %>%
  mutate(PC1 = ifelse(!is.na(Bars),round(transform_PCs(PC1),2),NA),
         PC2 = ifelse(!is.na(Bars),round(transform_PCs(PC2),2),NA),
         PC_d1 = ifelse(!is.na(Bars),round(transform_PCs(PC_d1),2),NA),
         ID = id) %>%
 select(-id) %>%
 select(ID,label:uni) %>%
 mutate(blue2 = ind+may)

write_delim(data_export,path = 'phenotypes.txt', delim = '\t')
