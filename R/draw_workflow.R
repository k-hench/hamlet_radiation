library(tidyverse)
library(ggraph)
library(ggiraph)
library(tidygraph)
library(cowplot)

clr <- RColorBrewer::brewer.pal(5, 'Set1')[c(3, 1)] %>% set_names(nm = c("start", 'product'))

read_dot <- function(file,
                     point_types = tibble(label = 'dummy', type = 'dummy')){
  raw_lines <- read_lines(file) %>%
    tibble(line = .) %>%
    filter(grepl('^p', line)) %>%
    mutate(line = line %>% 
             str_remove_all(pattern = '^ | $|;$|\\[.*label=\\"') %>% 
             str_remove(pattern = '\\"\\]$') %>%
             str_remove(pattern = ' \\[.*\\]'))
  
  nodes <- raw_lines %>%
    filter(!grepl('->', line)) %>%
    separate(line, into = c('point', 'label'), sep = ' ') %>%
    mutate(label = ifelse(is.na(label),'',label)) %>%
    group_by(point, label) %>%
    count() %>%
    ungroup() %>%
    left_join(., point_types, by = "label") %>%
    mutate(type = type %>% ifelse(point == 'p0',"start",.))
  
  edges <- raw_lines %>%
    filter(grepl('->', line)) %>%
    mutate(line = line %>% str_replace(pattern = " -> ", replacement = " ")) %>%
    separate(line, into = c('from', 'to', 'label'), sep = ' ') %>%
    mutate(label = ifelse(is.na(label),'',label)) %>%
    left_join(., point_types, by = "label")
  
  list(nodes = nodes, edges = edges)
}

dot_plot <- function(file, git,  point_types = tibble(label = 'dummy', type = 'dummy'),...){
  dot_data <- read_dot(file = file, 
                       point_types = point_types)
  
  g <- tbl_graph(nodes = dot_data$nodes, edges = dot_data$edges) %>%
    ggraph(layout = 'kk')+
    coord_flip()+
    annotation_custom(grid::textGrob(label = git,
                                     gp = grid::gpar(fontsize = 300,
                                                     fontface = 'bold',
                                                     col = '#E0e0e0')))+
    geom_edge_link(aes(label = label,color = type), 
                   angle_calc = 'along',
                   label_dodge = unit(2.5, 'mm'), 
                   start_cap = circle(8, 'pt'),
                   end_cap = circle(8, 'pt'),
                   arrow = arrow(type = 'closed',
                                 length = unit(2,'pt')))+
    geom_node_point(aes(color = type), shape = 1, size = 4) +
    scale_edge_color_manual(values = clr, na.value = 'black')+
    scale_color_manual(values = clr, na.value = 'black')+
    theme_void()+
    theme(legend.position = 'none')
  
  my_gg <- g + geom_point_interactive(aes(tooltip = label,x=x,y=y,color = type), size = 2) 
}

tibbler <- function(products){ tibble(label = products, type = 'product') }