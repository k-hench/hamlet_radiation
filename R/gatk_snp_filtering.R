library(tidyverse)

data <- read_delim('1_genotyping/1_raw_vcfs/raw_var_sites.vcf.gz.table.txt.gz', delim='\t') %>% mutate(logFS = log10(FS))

names(data)

vars <- c("MQ","QD","logFS","MQRankSum","ReadPosRankSum")
basedir <- 'figures/snp_filtering/'

rngs <- list(MQ = c(0,90),
             QD = c(0,45),
             logFS = c(-1,3),
             MQRankSum = c(-5,5),
             ReadPosRankSum = c(-5,5))

thresholds <- list(MQ = c(52,65),
                   QD = c(2.5),
                   logFS = c(log10(25)),
                   MQRankSum = c(-.2,.2),
                   ReadPosRankSum = c(-2,2))

plotting <- function( data, var, rng, thres ){
  p1 <- ggplot( data, aes_string( x = var )) +
    geom_density( fill = rgb(0, 0, .5, .5), n = 4096 ) +
    geom_vline( xintercept = thres )+
    coord_cartesian(xlim = rng) +
    scale_x_continuous(name = str_c(var,' ',str_c(thres, collapse = ', ')))
}

plts <- purrr::pmap(.l = tibble(var = vars,rng = rngs, thres = thresholds),
                      .f = plotting,data=data)
p1 <- cowplot::plot_grid(plotlist = plts,ncol = 1)

A4 <- c(210, 297)
ggsave(plot = p1, filename = 'figures/snp_filtering/gatk_filter_stats.pdf',
       width = .95 * A4[1], height = .95 * A4[2], 
       units ='mm', device = cairo_pdf)
