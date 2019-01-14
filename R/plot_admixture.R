#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=FALSE)
#args = list('admixture_report.bel.txt','pop.fake.txt','test.biallelic.phased','~/Desktop/chapter2/R/project_config.R', '.bel',8)
#args = list('admixture_report.bel.txt','pop.bel.txt','~/Desktop/chapter2/R/project_config.R', '.bel',8)
args = args[7:11]
print(args)
# config ----------------------
admx_rep <- as.character(args[1])
sample_id <- as.character(args[2])
config_script <- as.character(args[3])
loc <- as.character(args[4])
n <- as.numeric(args[5])

loc_clr <- 'black'

whisker <- .01
ybase <- 1.05  
ybase_loc <- 1.1
ylab_loc <- .02
# ----------------------
library(tidyverse)
library(hypoimg)
source(config_script)
fll <- fll_n(n)

import_admix <- function(K,loc){
  tbl <- read.table(str_c("hapmap",loc,'.',K,".Q",sep = '')) %>% 
    as_tibble() %>%
    set_names(.,nm = c(str_c('pop_',1:K))) %>% 
    mutate(ind_nr = row_number()) %>% 
    gather(key = 'pop',
           value = 'Ancestry',-ind_nr) %>%
    arrange(ind_nr) %>% 
    mutate(id = id_labs$id[ind_nr]) %>% 
    left_join(.,id_labs)
  #,
  #         id = factor(id,levels = id_labs$id[order(id_labs$order)])) %>%
  #  select(ind_nr,id,pop,Ancestry) 
  }

plot_admix <- function(K,loc){
  p <- ggplot(import_admix(K,loc), aes(x = order,y = Ancestry,fill = pop))
  
  if(n_locs > 1){
    p <- p + 
      geom_segment(inherit.aes = FALSE,
                   data = loc_box, aes(x = start, xend = end,
                                       y = ybase_loc, yend = ybase_loc),
                   col = loc_clr) +
      geom_segment(inherit.aes = FALSE,
                   data = loc_box, aes(x = end, xend = end,
                                       y = -whisker+ybase_loc,yend=whisker+ybase_loc),
                   col = loc_clr) +
      geom_segment(inherit.aes = FALSE, 
                   data = loc_box, aes(x = start, xend = start, 
                                       y = -whisker+ybase_loc, yend = whisker+ybase_loc),
                   col = loc_clr) +
      geom_label(inherit.aes = FALSE, 
                data = loc_box, aes(x = (start+end)/2, 
                                    y = ybase_loc+ylab_loc,label = loc_names[loc]),
                col = loc_clr, size = 2,label.padding = unit(2,'pt'),nudge_y = -.012) +
      scale_y_continuous(str_c("Ancestry, K: ",K),
                         expand = c(0.005,0.005), limits = c(0,1.15))
  } else {
    p <- p + scale_y_continuous(str_c("Ancestry, K: ",K), expand = c(0.005,0.005))
  }
  
  p <- p +
    geom_bar(stat = 'identity',col='black')+
    # ggtitle()+
    #  scale_fill_manual(values = scico::scico(3,palette = "lajolla"))+
    geom_segment(inherit.aes = FALSE,data = spec_box,aes(x=start,xend=end,y=ybase,yend=ybase,col=spec))+
    geom_segment(inherit.aes = FALSE,data = spec_box,aes(x=end,xend=end,y=-whisker+ybase,yend=whisker+ybase,col=spec))+
    geom_segment(inherit.aes = FALSE,data = spec_box,aes(x=start,xend=start,y=-whisker+ybase,yend=whisker+ybase,col=spec))+
    scale_color_manual(values = clr2, guide = FALSE)+
    scale_fill_manual(name = str_c('K: ',K,sep=''), values = fll)+
    scale_x_discrete(expand = c(0.005,0.005))+
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
          #element_text(angle=90)
          )
  
  p
  
}

plot_admix_pure <- function(K,loc){plot_admix(K,loc)+theme(legend.position = 'none')} 

id_labs <- read_delim(sample_id,delim = '\t',col_names = c('id','spec','loc')) %>% 
  mutate(order = str_c(loc,spec,id))

n_locs <- id_labs$loc %>% as.factor() %>% levels() %>% length()

loc_box <- id_labs %>%
  group_by(loc) %>% 
  count() %>% 
  ungroup() %>%
  mutate(start = cumsum(lag(x = n,default = 0))+.6,
         end = start+n-.2)

spec_box <- id_labs %>%
  mutate(order = str_c(loc,spec)) %>%
  group_by(order) %>% 
  count() %>% 
  ungroup() %>%
  mutate(start = cumsum(lag(x = n,default = 0))+.6,
         end = start+n-.2,
        spec = str_extract(string = order,'...$'))

data <- read_delim(admx_rep,delim=' ',
                   col_names = c('cv', 'error', 'k', 'value')) %>% 
  mutate(K = row_number()) %>% 
  filter(K<=n)

preplot <- ggplot(data, aes(x = K, y = value)) +
  geom_area(fill = rgb(.1, .1, .1, .1)) +
  geom_point(aes(fill = factor(K)),shape = 21, size=1.5) +
  scale_fill_manual(values = fll_fun(n), guide = FALSE) +
  theme_minimal()

sp_list <- spec_box$spec
legend_grob <- hypo_legend_single(species = sp_names[sp_list],
                                  color_map = clr2[sp_list],
                                  circle_color = 'black',
                                  plot_names = TRUE,
                                  circle_lwd = .5,
                                  ncol = min(length(sp_list),6)) %>%
  ggplotGrob()

legends <- plot_grid((plot_admix(n,loc = loc) + 
                        guides(fill = guide_legend(title = "Assigned Population:",nrow = 1)) +
                        theme(legend.position = 'bottom')) %>% 
                       get_legend(),
                     legend_grob,
                     ncol=1,
                     rel_heights = c(.4,1))

plot_list <- purrr::map(2:n, plot_admix_pure, loc = loc)
plots_composition <- cowplot::plot_grid(plotlist = rlist::list.prepend(plot_list, preplot), ncol = 2)
plot_complete <- cowplot::plot_grid(plots_composition,legends,ncol = 1,rel_heights = c(1,.15))

ggsave(plot_complete,filename = str_c('admixture',loc,'.pdf'), width = 8,height = 9)
