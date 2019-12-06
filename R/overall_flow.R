library(tidyverse)
library(ggraph)
library(tidygraph)
library(cowplot)

node_meta <- tibble(name = 0:9 %>% as.character(),
                   label = ifelse(name == 0, 'raw data',str_c('git ', name, '.x')) %>% 
                     str_c(c('','\ngenotyping','\ngenotyping all bp', '\nfst/gxp','\ndxy','\nphylo','\nrecombination','\nheterozygosity','\nmsmc','\nfigs')))

start_buffer <- tibble(from = 1:9, start_buffer = c(10, 28, rep(25,7)))
end_buffer <- tibble(from = 1:9, end_buffer = c(20, 25, rep(25,7)))


graph <- as_tbl_graph(
  tibble(
    from = c(0,1,1,2,1,1,1,1,3:8,3),
    to = c(1,2:8,rep(9,6),4),
    weight = 1
  )
)

graph <- graph %>%
  activate('nodes') %>%
  left_join(node_meta) %>%
  activate('edges') %>%
  left_join(start_buffer)%>%
  left_join(end_buffer) %>%
  left_join(node_meta %>% 
              mutate(from = name %>% as.integer()) %>% 
              select(from, label) %>% set_names(nm = c('from', 'from_label'))) %>%
  left_join(node_meta %>% 
            mutate(to = name %>% as.integer()) %>% 
            select(to, label) %>% set_names(nm = c('to', 'to_label'))) 

p_overall <- ggraph(graph,layout = 'sugiyama') + 
  coord_equal()+
  geom_edge_link(aes(start_cap = label_rect(from_label),
                     end_cap = circle(end_buffer, 'pt')),
                 alpha =.4,
                 angle_calc = 'along',
                 label_dodge = unit(2.5, 'mm'),
                 arrow = arrow(type = 'closed',
                               length = unit(3,'pt')))+
  geom_node_label(aes(label = label),
                  label.size = 0,fill = rgb(1,1,1,0))+
  scale_x_continuous(expand = c(.2,.3))+
  scale_y_continuous(expand = c(.05,.1))+
  theme_graph()

