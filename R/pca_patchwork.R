library(GenomicOriginsScripts)
library(hypogen)
library(patchwork)

types <- c("whg", "subset_non_diverged")

pca_dir <- "2_analysis/pca/"

titles <- c(loc_names, all = "All Samples", hamlets_only = "Hamlets only")
main_title <- c(whg = "Whole Genome (all SNPs)",
                subset_non_diverged = "SNP subset (diverged regions excluded)")

patchwork_pcas <- function(type){
  score_files <- dir(pca_dir, pattern = str_c(type,".scores"))
  exp_files <- dir(pca_dir, pattern = str_c(type,".exp_var"))
  exp_extension <- ".exp_var.txt.gz"
  
  read_pca <- function(scores){
    score_base <- scores %>% str_remove(pattern = ".scores.txt.gz")
    who <- ifelse(score_base == type, "all", str_remove(score_base, "\\..*$"))
    what <- score_base %>% str_remove("^.*\\.")
    
    score_tibble <- read_tsv(str_c(pca_dir, scores))
    exp_tibble <- read_tsv(str_c(pca_dir, score_base, exp_extension)) 
    tibble(who = who, what = what,
           scores = list(score_tibble), exp_var = list(exp_tibble) )
  }
  
  score_tibble <- score_files %>% map_dfr(read_pca) %>% 
    unnest(cols = c(scores))
  exp_tibble <- score_files %>% map_dfr(read_pca) %>% 
    unnest(cols = c(exp_var)) %>%
    select(-EV_nr) %>%
    pivot_wider(names_from = EV, values_from = exp_var)
  
  plot_pca <- function(who_plot){
    score_tibble %>%
      filter(who == who_plot) %>%
      mutate(spec = str_sub(id, -6, -4)) %>%
      ggplot(aes(x = EV01, y = EV02, fill = spec, color = spec))+
      geom_point( shape = 21, aes(color = after_scale(prismatic::clr_darken(fill)))) +
      ggforce::geom_mark_ellipse(expand = unit(5, "pt"), fill = "transparent", linetype = 3)+
      labs(title = titles[who_plot],
           x = str_c("EV01 (", sprintf("%.2f", exp_tibble$EV01[exp_tibble$who == who_plot]),"%)"),
           y = str_c("EV02 (", sprintf("%.2f", exp_tibble$EV02[exp_tibble$who == who_plot]),"%)"))+
      scale_fill_manual(values = clr) +
      scale_color_manual(values = clr2)
  }
  
  (plot_pca("bel") + plot_pca("hon") + plot_pca("pan"))/
    (plot_pca("all") | plot_pca("hamlets_only")) +
    plot_layout(heights = c(1/3, .5)) +
    plot_annotation(title = main_title[type]) &
    theme_minimal() &
    theme(text = element_text(size = 8),
          legend.position = "none",
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill = "transparent",
                                          color = "lightgray",
                                          size = .5))
}

p1 <- patchwork_pcas(types[[1]])
p2 <- patchwork_pcas(types[[2]])

cowplot::plot_grid(p1, p2)

# Searching floflo 
# clr_swap <- clr
# clr_swap[] <- "gray"
# clr_swap["flo"] <- "red"
# p1 & scale_color_manual(values = clr_swap)& scale_fill_manual(values = clr_swap)

scl <- .7
ggsave("~/Desktop/whg_pcas.pdf", width = 16 * scl, height = 9 * scl, device = cairo_pdf)
