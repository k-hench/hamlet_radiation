---
output: html_document
editor_options:
  chunk_output_type: console
---
# Supplementary Figure 2






## Summary

This is the accessory documentation of Figure S2.
The Figure can be recreated by running the **R** script `plot_SF2.R`:

```sh
cd $BASE_DIR

Rscript --vanilla R/fig/plot_SF2.R \
    ressources/Rabosky_etal_2018/
```

## Details of `plot_SF2.R`

In the following, the individual steps of the R script are documented.
It is an executable R script that depends on the accessory R package [**GenomicOriginsScripts**](https://k-hench.github.io/GenomicOriginsScripts), as well as on the packages [**hypoimg**](https://k-hench.github.io/hypoimg), [**hypogen**](https://k-hench.github.io/hypogen) and [**patchwork**](https://patchwork.data-imaginist.com/)

### Config

The scripts start with a header that contains copy & paste templates to execute or debug the script:


```r
#!/usr/bin/env Rscript
# run from terminal:
# Rscript --vanilla R/fig/plot_SF2.R \
#     ressources/Rabosky_etal_2018/
# ===============================================================
# This script produces Suppl. Figure 2 of the study "Rapid radiation in a
# highly diverse marine environment" by Hench, Helmkampf, McMillan and Puebla
# ---------------------------------------------------------------
# ===============================================================
# args <- c( "ressources/Rabosky_etal_2018/" )
# script_name <- "R/fig/plot_SF2.R"
args <- commandArgs(trailingOnly = FALSE)
```

The next section processes the input from the command line.
It stores the arguments in the vector `args`.
The needed R packages are loaded and the script name and the current working directory are stored inside variables (`script_name`, `plot_comment`).
This information will later be written into the meta data of the figure to help us tracing back the scripts that created the figures in the future.

Then we drop all the imported information besides the arguments following the script name and print the information to the terminal.


```r
# setup -----------------------
renv::activate()
library(BAMMtools)
library(GenomicOriginsScripts)
library(ggplotify)
library(patchwork)
library(ggforce)
library(glue)
library(ggtext)
library(hypoimg)

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
```

```r
#> ── Script: R/fig/plot_SF2.R ────────────────────────────────────────────
#> Parameters read:
#> ★ 1: ressources/Rabosky_etal_2018/
#> ────────────────────────────────────────── /current/working/directory ──
```

The directory containing the PCA data is received and stored in a variable.
Also the default color scheme is updated and the size of the hamlet ann.


```r
# config -----------------------
basepath <-  as.character(args[1])

### Evolutionary rate analysis of Serraninae
### ----------------------------------------
### Code adapted from bamm-project.org and Rabosky et al. 2018 (Nature)

source(paste(basepath, "scripts/supporting_fxns/PlottingFunctions.R", sep = ""))
```



```r
## Import FToL data
eventfile_vr <- paste(basepath, "dataFiles/bamm_results/12k_tv1/event_data_thinned.csv", sep="")
treefile <- paste(basepath, "dataFiles/bamm_results/12k_tv1/bigfish_no_outgroup.tre", sep="")
```



```r
tree_ftol <- read.tree(treefile)

## Map event data onto time-calibrated tree
bamm_ftol <- getEventData(tree_ftol, eventfile_vr, burnin = 0)
```



```r
## Extract Serraninae subtree
Serraninae <- c("Hypoplectrus_gemma", "Hypoplectrus_unicolor", "Hypoplectrus_gummigutta", "Hypoplectrus_chlorurus", "Hypoplectrus_aberrans", "Hypoplectrus_nigricans",
                "Hypoplectrus_guttavarius", "Hypoplectrus_indigo", "Hypoplectrus_puella", "Serranus_tortugarum", "Serranus_tabacarius", "Schultzea_beta",
                "Diplectrum_formosum", "Diplectrum_bivittatum", "Diplectrum_pacificum", "Diplectrum_maximum", "Serranus_notospilus", "Serranus_phoebe",
                "Serranus_psittacinus", "Serranus_baldwini", "Serranus_tigrinus", "Paralabrax_albomaculatus", "Paralabrax_dewegeri", "Paralabrax_callaensis",
                "Paralabrax_loro", "Paralabrax_auroguttatus", "Paralabrax_clathratus", "Paralabrax_humeralis", "Paralabrax_nebulifer", "Paralabrax_maculatofasciatus",
                "Zalanthias_kelloggi", "Serranus_cabrilla", "Serranus_atricauda", "Serranus_scriba", "Serranus_hepatus", "Serranus_accraensis", "Centropristis_striata",
                "Chelidoperca_occipitalis", "Chelidoperca_investigatoris", "Chelidoperca_pleurospilus")

Hamlets <- c("Hypoplectrus_gemma", "Hypoplectrus_unicolor", "Hypoplectrus_gummigutta", "Hypoplectrus_chlorurus", "Hypoplectrus_aberrans", "Hypoplectrus_nigricans",
             "Hypoplectrus_guttavarius", "Hypoplectrus_indigo", "Hypoplectrus_puella")
```



```r
bamm_serrn <- subtreeBAMM(bamm_ftol, tips = Serraninae)
tree_serrn <- as.phylo(bamm_serrn)
```



```r
## Mean phylorate plot
bamm_serrn_abbr <- bamm_serrn
bamm_serrn_abbr$tip.label <- bamm_serrn$tip.label %>%
  str_replace(pattern = "([A-Z])[a-z]*_([a-z]*)", "italic(\\1.~\\2)") %>%
  str_replace(pattern = "C.", "Ch.") %>%
  str_replace(pattern = "Ch.~striata", "Cp.~striata") %>%
  str_replace(pattern = "S.", "Se.") %>%
  str_replace(pattern = "Se.~beta", "Sc.~'beta'") %>%
  str_replace(pattern = "P.", "Pa.") %>%
  str_replace(pattern = "Z.", "Pl.") %>%
  ggplot2:::parse_safe()
```



```r
## Credible sets of shift configurations
css <- credibleShiftSet(bamm_serrn_abbr, expectedNumberOfShifts = 1, threshold = 5, set.limit = 0.95)
summary(css)
```



```r
clr_tree <- scico::scico(6, palette = "berlin") %>%
  prismatic::clr_desaturate(shift = .4) %>%
  prismatic::clr_darken(shift = .2)

clr_tree2 <- colorRampPalette(RColorBrewer::brewer.pal(9,"RdYlBu"))(64) %>% rev()

clr_shift <- "red"

css2 <- css
css2$marg.probs["47"] <- .1
```



```r
p1 <- as.grob(function(){
  par(mar = c(0,0,0,0))
  plot.credibleshiftset(css2, logcolor = TRUE,
                        add.freq.text = FALSE,
                        border = FALSE,
                        shiftColor = clr_shift,
                        lwd = 2,
                        labels = TRUE,
                        legend = FALSE,
                        pal = clr_tree, cex = .5)

  leg_shift_x <- 1.3
  leg_shift_y <- 5
  text(x = c(21.2, 40.2), y = c(15.6, 33.25),
       label = "\U2605", family = "DejaVu Sans", col = clr_shift, cex = .5)
  lines(x = c(0,25) + leg_shift_x,
        y = c(1.5, 1.5) + leg_shift_y,
        col = "darkgray")
  text(x = 12.5 + leg_shift_x,
       y = .5 + leg_shift_y,
       labels = "25 MYR",
       cex = .4,
       col = "darkgray")
})
```



```r
## Macroevolutionary cohort analysis
cmat <- getCohortMatrix(bamm_serrn)
p2 <- as.grob(function(){
cohorts(cmat, bamm_serrn,
        lwd = 1.5,
        labels = FALSE,
        legend = FALSE,
        ofs = 0,
        use.plot.bammdata = TRUE,
        pal = clr_tree,
        col = clr_tree2,
        cex.axis = 0.1)
})
```



```r
p_done <- (ggplot() +
              geom_point(data = tibble(v = c(.056, 2.4)),
                         x = .5, y = .5, aes(color = v),alpha = 0) +
              scale_color_gradientn("Speciation Rate",
                                    colours = clr_tree,
                                    limits = c(.056, 2.4)) +
              geom_richtext(data = tibble(x = -.05,
                                          y = .12,
                                          lab = glue("Posterior Frequency: {css$frequency}<br>Marginal Shift Prob.: {css$marg.probs['47']}")),
                            aes(x = x, y = y, label = lab),
                            size = plot_text_size_small / .pt,
                            color = clr_shift,
                            hjust = 0,
                            label.size = 0,
                            label.color = "transparent")+
              geom_bezier0(data = tibble(x = c(.62,.5,.32), y = c(.185,.12,.12)),
                           aes(x,y, group = 1),
                           size = .3,
                           color = prismatic::clr_alpha(clr_shift,.3))+
              annotation_custom(p1,xmin = -.3, ymin = -.2,
                                xmax = 1.1, ymax = 1.13)  +
    (ggplot() +
       geom_point(data = tibble(v = c(0, 1)),
                  x = .5, y = .5, aes(color = v),alpha = 0)+
       scale_color_gradientn("Pairwise Correlation", colours = clr_tree2, limits = c(0, 1))+
       annotation_custom(p2, xmin = -.2, ymin = -.25,
                         xmax = 1.2, ymax = 1.075) )  &
    plot_annotation(tag_levels = "a") &
    coord_cartesian(xlim = c(0,1),
                    ylim = c(0,1)) &
      guides(color = guide_colorbar(title.position = "top",
                                    direction = "horizontal",
                                    barheight = unit(3, "pt"),
                                    barwidth = unit(100, "pt"),
                                    ticks.colour = "white"))) &
    theme_minimal(base_size = plot_text_size) &
    theme(legend.position = c(.5, -.03),
          legend.justification = c(.5, 0),
          legend.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.tag = element_text(hjust = 0),
          panel.background = element_blank(),
          plot.background = element_blank())
```

Finally, we can export Figure S2.


```r
hypo_save("figures/SF2.pdf",
       plot = p_done,
       width = f_width,
       height = f_width*.5,
       comment = plot_comment,
       device = cairo_pdf,
       bg = "transparent")
```
