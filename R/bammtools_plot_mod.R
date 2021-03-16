bammplot_k <- function (x, tau = 0.01, method = "phylogram", xlim = NULL, 
                        ylim = NULL, vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, 
                        legend = FALSE, spex = "s", lwd = 1, cex = 1, pal = "RdYlBu", 
                        mask = integer(0), mask.color = gray(0.5), colorbreaks = NULL, 
                        logcolor = FALSE, breaksmethod = "linear", color.interval = NULL, 
                        JenksSubset = 20000, par.reset = FALSE, direction = "rightwards",
                        labelcolor = "black", 
                        ...) 
{
  if (inherits(x, "bammdata")) {
    if (attributes(x)$order != "cladewise") {
      stop("Function requires tree in 'cladewise' order")
    }
    phy <- BAMMtools:::as.phylo.bammdata(x)
  }
  else stop("Object ephy must be of class bammdata")
  if (!spex %in% c("s", "e", "netdiv")) {
    stop("spex must be 's', 'e' or 'netdiv'.")
  }
  if (length(pal) == 1 && !pal %in% names(get("palettes", 
                                              envir = .colorEnv)) && pal != "temperature" && pal != 
      "terrain") 
    pal <- rep(pal, 3)
  else if (length(pal) == 2) 
    pal <- c(pal, pal[2])
  if (breaksmethod == "linear" & !is.null(color.interval)) {
    if (length(color.interval) != 2) {
      stop("color.interval must be a vector of 2 numeric values.")
    }
  }
  if (!is.binary.phylo(phy)) {
    stop("Function requires fully bifurcating tree")
  }
  if (any(phy$edge.length == 0)) {
    warning("Tree contains zero length branches. Rates for these will be NA and coerced to zero")
  }
  if (!("dtrates" %in% names(x))) {
    x <- dtRates(x, tau)
  }
  NCOLORS <- 64
  if (!is.null(color.interval)) {
    if (x$type == "trait") {
      ratesRange <- range(x$dtrates$rates)
    }
    else if (x$type == "diversification") {
      if (tolower(spex) == "s") {
        ratesRange <- range(x$dtrates$rates[[1]])
      }
      else if (tolower(spex) == "e") {
        ratesRange <- range(x$dtrates$rates[[2]])
      }
      else if (tolower(spex) == "netdiv") {
        ratesRange <- range(x$dtrates$rates[[1]] - x$dtrates$rates[[2]])
      }
    }
    if (all(!is.na(color.interval))) {
      brks <- seq(min(color.interval[1], ratesRange[1]), 
                  max(color.interval[2], ratesRange[2]), length.out = (NCOLORS + 
                                                                         1))
      intervalLength <- length(which.min(abs(color.interval[1] - 
                                               brks)):which.min(abs(color.interval[2] - brks)))
    }
    else if (is.na(color.interval[1])) {
      brks <- seq(ratesRange[1], max(color.interval[2], 
                                     ratesRange[2]), length.out = (NCOLORS + 1))
      intervalLength <- length(1:which.min(abs(color.interval[2] - 
                                                 brks)))
    }
    else if (is.na(color.interval[2])) {
      brks <- seq(min(color.interval[1], ratesRange[1]), 
                  ratesRange[2], length.out = (NCOLORS + 1))
      intervalLength <- length(which.min(abs(color.interval[1] - 
                                               brks)):length(brks))
    }
    NCOLORS <- round((NCOLORS^2)/intervalLength)
  }
  if (is.null(colorbreaks)) {
    colorbreaks <- assignColorBreaks(x$dtrates$rates, NCOLORS, 
                                     spex, logcolor, breaksmethod, JenksSubset)
  }
  if (x$type == "trait") {
    colorobj <- BAMMtools:::colorMap(x$dtrates$rates, pal, colorbreaks, 
                                     logcolor, color.interval)
  }
  else if (x$type == "diversification") {
    if (tolower(spex) == "s") {
      colorobj <- BAMMtools:::colorMap(x$dtrates$rates[[1]], pal, 
                                       colorbreaks, logcolor, color.interval)
    }
    else if (tolower(spex) == "e") {
      colorobj <- BAMMtools:::colorMap(x$dtrates$rates[[2]], pal, 
                                       colorbreaks, logcolor, color.interval)
    }
    else if (tolower(spex) == "netdiv") {
      colorobj <- BAMMtools:::colorMap(x$dtrates$rates[[1]] - x$dtrates$rates[[2]], 
                                       pal, colorbreaks, logcolor, color.interval)
    }
  }
  else {
    stop("Unrecognized/corrupt bammdata class. Type does not equal 'trait' or 'diversification'")
  }
  edge.color <- colorobj$cols
  tH <- max(x$end)
  phy$begin <- x$begin
  phy$end <- x$end
  tau <- x$dtrates$tau
  if (method == "polar") {
    ret <- setPolarTreeCoords(phy, vtheta, rbf)
    rb <- tH * rbf
    p <- BAMMtools:::mkdtsegsPolar(ret$segs[-1, ], tau, x$edge)
  }
  else if (method == "phylogram") {
    ret <- BAMMtools:::setPhyloTreeCoords(phy)
    p <- BAMMtools:::mkdtsegsPhylo(ret$segs[-1, ], tau, x$edge)
  }
  else {
    stop("Unimplemented method")
  }
  x0 <- c(ret$segs[1, 1], p[, 1])
  x1 <- c(ret$segs[1, 3], p[, 2])
  y0 <- c(ret$segs[1, 2], p[, 3])
  y1 <- c(ret$segs[1, 4], p[, 4])
  offset <- table(p[, 5])[as.character(unique(p[, 5]))]
  if (length(mask)) {
    edge.color[p[, 5] %in% mask] <- mask.color
  }
  arc.color <- c(edge.color[1], edge.color[match(unique(p[, 
                                                          5]), p[, 5]) + offset])
  edge.color <- c(edge.color[1], edge.color)
  if (show) {
    op <- par(no.readonly = TRUE)
    if (length(list(...))) {
      par(...)
    }
    if (legend) {
      par(mar = c(5, 4, 4, 5))
    }
    plot.new()
    ofs <- 0
    if (labels) {
      if (method == "phylogram") 
        ofs <- max(nchar(phy$tip.label) * 0.03 * cex * 
                     tH)
      else ofs <- max(nchar(phy$tip.label) * 0.03 * cex)
    }
    if (method == "polar") {
      if (is.null(xlim) || is.null(ylim)) {
        if (is.null(xlim)) 
          xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs)
        if (is.null(ylim)) 
          ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs)
      }
      plot.window(xlim = xlim, ylim = ylim, asp = 1)
      segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, 
               lend = 2)
      arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + 
                                                  phy$end/tH), border = arc.color, lwd = lwd)
      if (labels) {
        for (k in 1:length(phy$tip.label)) {
          text(ret$segs[-1, ][phy$edge[, 2] == k, 3], 
               ret$segs[-1, ][phy$edge[, 2] == k, 4], phy$tip.label[k], 
               cex = cex, srt = (180/pi) * ret$arcs[-1, 
               ][phy$edge[, 2] == k, 1], adj = c(0, NA),
               col = labelcolor[[k]])
        }
      }
    }
    if (method == "phylogram") {
      direction <- match.arg(direction, c("rightwards", 
                                          "leftwards", "downwards", "upwards"))
      if (direction == "rightwards") {
        bars <- BAMMtools:::redirect(cbind(x0, y0, x1, y1), 0)
        arcs <- BAMMtools:::redirect(ret$arcs, 0)
        bars[, c(1, 3)] <- tH * bars[, c(1, 3)]
        arcs[, c(1, 3)] <- tH * arcs[, c(1, 3)]
        ret$segs[-1, c(1, 3)] <- tH * ret$segs[-1, c(1, 
                                                     3)]
      }
      else if (direction == "leftwards") {
        bars <- BAMMtools:::redirect(cbind(x0, y0, x1, y1), pi)
        bars[, c(2, 4)] <- abs(bars[, c(2, 4)])
        arcs <- BAMMtools:::redirect(ret$arcs, pi)
        arcs[, c(2, 4)] <- abs(arcs[, c(2, 4)])
        bars[, c(1, 3)] <- tH * bars[, c(1, 3)]
        arcs[, c(1, 3)] <- tH * arcs[, c(1, 3)]
        ret$segs[-1, c(1, 3)] <- -tH * ret$segs[-1, 
                                                c(1, 3)]
      }
      else if (direction == "downwards") {
        bars <- BAMMtools:::redirect(cbind(x0, y0, x1, y1), -pi/2)
        arcs <- BAMMtools:::redirect(ret$arcs, -pi/2)
        bars[, c(2, 4)] <- tH * bars[, c(2, 4)]
        arcs[, c(2, 4)] <- tH * arcs[, c(2, 4)]
        ret$segs <- BAMMtools:::redirect(ret$segs, -pi/2)
        ret$segs[, c(2, 4)] <- tH * ret$segs[, c(2, 
                                                 4)]
      }
      else if (direction == "upwards") {
        bars <- BAMMtools:::redirect(cbind(x0, y0, x1, y1), pi/2)
        bars[, c(1, 3)] <- abs(bars[, c(1, 3)])
        arcs <- BAMMtools:::redirect(ret$arcs, pi/2)
        arcs[, c(1, 3)] <- abs(arcs[, c(1, 3)])
        bars[, c(2, 4)] <- tH * bars[, c(2, 4)]
        arcs[, c(2, 4)] <- tH * arcs[, c(2, 4)]
        ret$segs <- BAMMtools:::redirect(ret$segs, pi/2)
        ret$segs[, c(1, 3)] <- abs(ret$segs[, c(1, 3)])
        ret$segs[, c(2, 4)] <- tH * ret$segs[, c(2, 
                                                 4)]
      }
      if (is.null(xlim) && direction == "rightwards") 
        xlim <- c(0, tH + ofs)
      if (is.null(xlim) && direction == "leftwards") 
        xlim <- c(-(tH + ofs), 0)
      if (is.null(ylim) && (direction == "rightwards" || 
                            direction == "leftwards")) 
        ylim <- c(0, phy$Nnode)
      if (is.null(xlim) && (direction == "upwards" || 
                            direction == "downwards")) 
        xlim <- c(0, phy$Nnode)
      if (is.null(ylim) && direction == "upwards") 
        ylim <- c(0, tH + ofs)
      if (is.null(ylim) && direction == "downwards") 
        ylim <- c(-(tH + ofs), 0)
      plot.window(xlim = xlim, ylim = ylim)
      segments(bars[-1, 1], bars[-1, 2], bars[-1, 3], 
               bars[-1, 4], col = edge.color[-1], lwd = lwd, 
               lend = 2)
      isTip <- phy$edge[, 2] <= phy$Nnode + 1
      isTip <- c(FALSE, isTip)
      segments(arcs[!isTip, 1], arcs[!isTip, 2], arcs[!isTip, 
                                                      3], arcs[!isTip, 4], col = arc.color[!isTip], 
               lwd = lwd, lend = 2)
      if (labels) {
        if (direction == "rightwards") 
          text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
               phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
               pos = 4, offset = 0.25,
               col = labelcolor)
        else if (direction == "leftwards") 
          text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
               phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
               pos = 2, offset = 0.25,
               col = labelcolor)
        else if (direction == "upwards") 
          text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
               phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
               pos = 4, srt = 90, offset = 0,
               col = labelcolor)
        else if (direction == "downwards") 
          text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
               phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
               pos = 2, srt = 90, offset = 0,
               col = labelcolor)
      }
    }
  }
  index <- order(as.numeric(rownames(ret$segs)))
  if (show) {
    if (method == "phylogram") {
      assign("last_plot.phylo", list(type = "phylogram", 
                                     direction = direction, Ntip = phy$Nnode + 1, 
                                     Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 
                                                                                       3], yy = ret$segs[index, 4], pp = par(no.readonly = TRUE)), 
             envir = .PlotPhyloEnv)
    }
    else if (method == "polar") {
      assign("last_plot.phylo", list(type = "fan", Ntip = phy$Nnode + 
                                       1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 
                                                                                            3], yy = ret$segs[index, 4], theta = ret$segs[index, 
                                                                                                                                          5], rb = rb, pp = par(no.readonly = TRUE)), 
             envir = .PlotPhyloEnv)
    }
    if (legend) {
      addBAMMlegend(x = list(coords = ret$segs[-1, ], 
                             colorbreaks = colorobj$breaks, palette = colorobj$colpalette, 
                             colordens = colorobj$colsdensity), location = "right")
    }
  }
  if (par.reset) {
    par(op)
  }
  invisible(list(coords = ret$segs[-1, ], colorbreaks = colorobj$breaks, 
                 palette = colorobj$colpalette, colordens = colorobj$colsdensity))
}
