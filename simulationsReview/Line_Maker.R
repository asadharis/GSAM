
# A function to create the default colors of ggplot.
#
# n: Number of different colors to generate. E.g. if you have
#   4 different colors then ggplot uses ggplot_cols(4).
ggplot_cols <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# This function creates a line with symbol which can be used in
# plot captions in Latex.
#
#
# pch: pch of symbol, no symbol drawn if NA.
# col: color passed on to `col` argument, can be string, integer or
#      color code as in output of ggplot_cols function.
# lty: type of line.
# col.label: (Optional) character can relabel color names to
#           be more informative. e.g. col.label = c("Experiment 1", "Experiment 2")
#
# Example: myf(NA, "red", 1) creates a solid red line.
#          Latex code to include figure:
#         \protect\includegraphics[scale=0.4]{NA_red_1.png}
#
myf <- function(pch = NA, col = "red", lty = 2,
                col.label = NA) {
  if(!is.na(col.label)) {
    name <- paste0(pch,"_", col.label, "_", lty)
  } else {
    name <- paste0(pch,"_", col, "_", lty)
  }

  png(paste0(name, ".png"), width = 1.4, height = 0.30, units = "in", res = 6000)
  par(mar = c(0,0,0,0))
  plot(0:1,c(0.4,0.6), type = "n", axes = FALSE, ann = FALSE)
  lines(0:1, rep(0.5,2), lwd = 3.5, col = col, lty = lty)
  if(!is.na(pch)) {
    points(c(0.33, 0.67), c(0.5,0.5), pch = pch, col = col, cex = 2.5)
  }
  dev.off()
}


cols <- ggplot_cols(2)
myf(pch = 16, col = cols[1], lty = 1, col.label = "coupled")
myf(pch = 17, col = cols[2], lty = 2, col.label = "decoupled")
