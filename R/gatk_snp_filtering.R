library(tidyverse)


data <- read_delim('1_genotyping/1_raw_vcfs/raw_var_sites.vcf.gz.table.txt.gz', delim='\t') %>% mutate(logFS = log10(FS))

names(data)

vars <- c("MQ","QD","logFS","MQRankSum","ReadPosRankSum")
basedir <- 'figures/snp_filtering/'

rngs <- list(c(0,90),
             c(0,45),
             c(-1,3),
             c(-5,5),
             c(-5,5))

plotting <- function(data,var,rng){
  p1 <- ggplot(data, aes_string( x = var )) +
    geom_density(fill=rgb(0,0,.5,.5)) +
    coord_cartesian(xlim = rng) 
}

plts <- purrr::map2(.x = vars,.y = rngs,.f = plotting,data=data)
p1 <- cowplot::plot_grid(plotlist = plts,ncol = 1)

A4 <- c(210, 297)
ggsave(plot = p1, filename = 'figures/snp_filtering/gatk_filter_stats.pdf',
       width = .95 * A4[1], height = .95 * A4[2], 
       units ='mm', device = cairo_pdf)
