clr <- c(
  abe = "#E5E5A1",
  flo = "#ABA7C4",
  gum = "#E3A258",
  ind = "#1C6FCE",
  may = "#7EA7C2",
  nig = "#000000",
  pue = "#E17366",
  ran = "#7EBDB3",
  tab = "#A1C75E",
  tor = "#E2B8CE",
  uni = "#ffffff"
)

clr2 <- c(
  abe = "#E5E5A1",
  flo = "#ABA7C4",
  gum = "#E3A258",
  ind = "#1C6FCE",
  may = "#7EA7C2",
  nig = "#000000",
  pue = "#E17366",
  ran = "#7EBDB3",
  tab = "#A1C75E",
  tor = "#E2B8CE",
  uni = "#CCCCCC"
)

plot_clr <- rgb(.2,.2,.2)
plot_size <- .2
outlr_clr <- rgb(1,0,0,.2)

loc_names <- c(
	bel = "Belize",
	hon = "Honduras",
	pan = "Panama"
)

clr_loc = c(
	bel = "#440154FF",
	hon = "#21908CFF",
	pan = "#FDE725FF"
)

fll_fun <- function(n){viridis::inferno(n)}
fll_n <- function(n){fll_fun(n) %>% setNames(., nm = str_c("pop_",1:n)) }

sp_labs <- c(
  abe = expression(italic(H.~aberrans)),
  flo = expression(italic(H.~floridae)),
  gum = expression(italic(H.~gummigutta)),
  ind = expression(italic(H.~indigo)),
  may = expression(italic(H.~maya)),
  nig = expression(italic(H.~nigricans)),
  pue = expression(italic(H.~puella)),
  ran = expression(italic(H.~randallorum)),
  tab = expression(italic(S.~tabacarius)),
  tor = expression(italic(S.~tortugarum)),
  uni = expression(italic(H.~unicolor))
)

sp_names <- c(
  abe = "aberrans",
  flo = "floridae",
  gum = "gummigutta",
  ind = "indigo",
  may = "maya",
  nig = "nigricans",
  pue = "puella",
  ran = "randallorum",
  tab = "tabacarius",
  tor = "tortugarum",
  uni = "unicolor"
)

shps <- c(bel = 21,
          hon = 23,
          pan = 22,
          flo = 24)

loc_labs <- c( bel = "Belize",
               hon = "Honduras",
               pan = "Panama")

project_case <- function(x){
	str_to_lower(x)
}

project_inv_case <- function(x){
  str_to_upper(x)
}