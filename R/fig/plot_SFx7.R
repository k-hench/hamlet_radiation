### Evolutionary rate analysis of Serraninae
### ----------------------------------------

### Code adapted from bamm-project.org and Rabosky et al. 2018 (Nature)
### Data from Rabosky et al. 2018 (Nature) -- Dryad: https://doi. org/10.5061/dryad.fc71cp4)
### Jul 30, 2021, Martin Helmkampf & Kosmas Hench

library(BAMMtools)
library(GenomicOriginsScripts)
library(ggplotify)
library(patchwork)

basepath <- "ressources/Rabosky_etal_2018/"

source(paste(basepath, "scripts/supporting_fxns/PlottingFunctions.R", sep=""))


## Import FToL data
eventfile_vr <- paste(basepath, "dataFiles/bamm_results/12k_tv1/event_data_thinned.csv", sep="")
treefile <- paste(basepath, "dataFiles/bamm_results/12k_tv1/bigfish_no_outgroup.tre", sep="")

tree_ftol <- read.tree(treefile)


## Map event data onto time-calibrated tree
bamm_ftol <- getEventData(tree_ftol, eventfile_vr, burnin = 0)


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

bamm_serrn <- subtreeBAMM(bamm_ftol, tips = Serraninae)
tree_serrn <- as.phylo(bamm_serrn)


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

# plot(bamm_serrn_abbr, method = "phylogram", spex = "s", breaksmethod = "linear", logcolor = TRUE, lwd = 3, labels = TRUE, legend = TRUE)


## Number of rate shifts
# summary(bamm_serrn)


## Credible sets of shift configurations
css <- credibleShiftSet(bamm_serrn_abbr, expectedNumberOfShifts = 1, threshold = 5, set.limit = 0.95)
summary(css)

clr_tree <- scico::scico(6, palette = "berlin") %>% 
  prismatic::clr_desaturate(shift = .4) %>% 
  prismatic::clr_darken(shift = .2)

clr_tree2 <- scales::colour_ramp(colors = RColorBrewer::brewer.pal(5,"PuOr")[c(5,4,3,2,1)])(seq(0,1,length.out = 501)) %>% 
  prismatic::clr_desaturate(shift = .2) %>% 
  prismatic::clr_lighten(shift = .2)

p1 <- as.grob(function(){
  par(mar = c(0,0,0,0))
  plot.credibleshiftset(css, logcolor = TRUE,
                        border = FALSE,
                        shiftColor = clr_tree2[[501]], lwd = 2, labels = TRUE, legend = FALSE,
                        pal = clr_tree, cex = .5)
})



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
        # col = clr_tree2, 
        cex.axis = .5)
})

p_done <- (ggplot() +
  annotation_custom(p1,xmin = -.3, ymin = -.33,xmax = 1.1,ymax = 1.04)) +
  (ggplot() +
     annotation_custom(p2,xmin = -.2, ymin = -.3))  &
  plot_annotation(tag_levels = "a") &
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) &
  theme(plot.tag = element_text(hjust = 0),
        panel.background = element_blank(),
        plot.background = element_blank())

ggsave("figures/SFx7.pdf", 
       plot = p_done, device = cairo_pdf,
       width = f_width, height = f_width*.5,
       bg = "transparent")
