clr <- c(
  abe = "#E5E5A1",
  flo = "#ABA7C4",
  gum = "#E3A258",
  ind = "#739FBD",
  may = "#C3C3C3",
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
  ind = "#739FBD",
  may = "#C3C3C3",
  nig = "#000000",
  pue = "#E17366",
  ran = "#7EBDB3",
  tab = "#A1C75E",
  tor = "#E2B8CE",
  uni = "#CCCCCC"
)

plot_clr <- rgb(.2,.2,.2)

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

v6 <- function(x){viridis::viridis(6)[x]}

fll = c( pop_1 = v6(1),
         pop_2 = v6(2),
         pop_3 = v6(3),
         pop_4 = v6(4),
         pop_5 = v6(5),
         pop_6 = v6(6) )

sp_labs <- c(
  abe = expression(italic(H.~aberrans)),
  flo = expression(italic(H.~floridae)),
  gum = expression(italic(H.~gumigutta)),
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
  gum = "gumigutta",
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