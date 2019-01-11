ceiling2d <- function(x){
  x <- x * 100
  x <- ceiling(x)
  x <- x / 100
  x
}

# function to import fst results (for use with purrr::pmap)
get_data <- function(file,pop1,pop2,base_dir){
  data_windows <- hypo_import_windows(file = str_c(base_dir,file), gz = TRUE) %>% 
    mutate(run = str_c(pop1,pop2,sep = '-'),
           window = 'bolditalic(F[ST])') %>% 
    group_by(CHROM) 
}

# helper function to directly get species comparison grobs 
# (for use with purrr::pmap)
plot_pair <- function(loc_run,loc,left,right,circle_fill_left,circle_fill_right){

  tibble(  loc_run = loc_run,loc = loc,
          grob = list(hypo_anno_pair_split(left,right,circle_fill_left,circle_fill_right,
                                           circle_color = 'black',circle_lwd = .2) %>% 
                        ggplotGrob()))
}

# scaling funtion to rescale global fst values for
# the horizontal bars on the side
#rescale_fst <- function(fst){
#  start <- hypogen::hypo_karyotype$GEND %>% last()*1.01
#  end <- 6.3e+08
#  fst_max <- max(globals$weighted)
#  
#  scales::rescale(fst,from = c(0,fst_max), to = c(start,end))
#}

rescale_fst <- function(fst){
  start <- 0#hypogen::hypo_karyotype$GEND %>% last()*1.01
  end <- 1
  fst_max <- max(globals$weighted)
  
  scales::rescale(fst,from = c(0,fst_max), to = c(start,end))
}

# helper function to crate a tibble for global fst values using
# th erescaling function above (prep for geom_rect)
fst_bar_row <- function(fst,loc_run){
  tibble(xmin = rescale_fst(0),
         xmax = rescale_fst(fst),
         xmin_org = 0,
         xmax_org = fst,
         ymin = 0,
         ymax= .15,loc_run = loc_run)
}

# helper function to refactor runs according to global fst values
refactor <- function(self,globals){
  factor(as.character(self$loc_run),
         levels = c(levels(globals$loc_run)))
}