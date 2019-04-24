sample_graph <- function(ID, n_sps = 50){
  cat(str_c('---',crayon::red(ID),'---\n'))

  system(str_c('gunzip -c ', vcf_file, ' | java -jar ', java_downsample, ' -n ',n_sps,' | gzip > downsample.vcf.gz'))

  popgr <- read.vcfR('downsample.vcf.gz') %>%
    vcfR2genind() %>%
    as.data.frame() %>%
    gstudio::to_mv() %>%
    popgraph::popgraph(x = ., groups = samples$grp)

  save(popgr,file = str_c('popgr.background.',ID,'.',n_sps,'.rda'))
  system(str_c('rm downsample.vcf.gz'))
}

generate_random_subset_graphs <- function(n_graphs, n_nsps){
  1:n_graphs %>%
    str_pad(width = 2,pad = '0') %>%
    purrr::map(sample_graph, n_sps = n_nsps)
}

slide_graph <- function(lg, start, end, n_snps, win_id, gwin_id){
  popgr <- read.vcfR(str_c('slide.', lg, '.', n_snps, '.', win_id, '.', gwin_id, '.vcf.gz')) %>%
    vcfR2genind() %>%
    as.data.frame() %>%
    gstudio::to_mv() %>%
    popgraph::popgraph(x = ., groups = samples$grp)

  save(popgr,file = str_c('popgr.', lg, '.', n_snps, '.', win_id, '.', gwin_id,'.rda'))
}

summarise_popgr <- function(file){
  graph_type <- str_split(file,pattern = '\\.') %>% unlist() %>% .[2]
  graph_idx <- str_split(file,pattern = '\\.') %>% unlist() %>% .[3]
  graph_snps <- str_split(file,pattern = '\\.') %>% unlist() %>% .[4]

  load(file)

  tibble(pop = V(popgr)$name,
         closeness = closeness(popgr),
         betweenness = betweenness(popgr),
         degree = degree( popgr ),
         strength = igraph::strength(popgr),
         eigenCent = evcent( popgr )$vector) %>%
    select(-pop) %>%
    summarise_all(mean) %>%
    mutate(n_nodes = igraph::vcount(popgr),
           n_edges = igraph::ecount(popgr),
           connected = igraph::is.connected(popgr),
           nr_cluster = igraph::no.clusters(popgr),
           radius = igraph::radius(popgr),
           mean_distance = igraph::mean_distance(popgr),
           sum_strength = sum(igraph::strength(popgr)),
           graph_type = graph_type,
           graph_idx = graph_idx,
           graph_snps = graph_snps) %>%
    select(graph_type:graph_snps,n_nodes:sum_strength,closeness:eigenCent)
}

summarise_slide_gr <- function(file){
  graph_lg <- str_split(file,pattern = '\\.') %>% unlist() %>% .[2]
  graph_snps <- str_split(file,pattern = '\\.') %>% unlist() %>% .[3]
  graph_win_id <- str_split(file,pattern = '\\.') %>% unlist() %>% .[4]
  graph_idx <- str_split(file,pattern = '\\.') %>% unlist() %>% .[5]

  load(file)

  tibble(pop = V(popgr)$name,
         closeness = closeness(popgr),
         betweenness = betweenness(popgr),
         degree = degree( popgr ),
         strength = igraph::strength(popgr),
         eigenCent = evcent( popgr )$vector) %>%
    select(-pop) %>%
    summarise_all(mean) %>%
    mutate(n_nodes = igraph::vcount(popgr),
           n_edges = igraph::ecount(popgr),
           connected = igraph::is.connected(popgr),
           nr_cluster = igraph::no.clusters(popgr),
           radius = igraph::radius(popgr),
           mean_distance = igraph::mean_distance(popgr),
           sum_strength = sum(igraph::strength(popgr)),
           lg = graph_lg,
           gwin_id = graph_idx,
           graph_win_id = graph_win_id,
           graph_snps = graph_snps) %>%
    select(lg:graph_snps,n_nodes:sum_strength,closeness:eigenCent)
}