gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- gg_color_hue(7)

plot(rep(1,7), 1:7, ylim = c(0.5,7.5), type = "n")
abline(h = 1, col = cols[1], lwd = 4, lty = 2)
abline(h = 2, col = cols[2], lwd = 4, lty = 2)
abline(h = 3, col = cols[3], lwd = 4, lty = 2)
abline(h = 4, col = cols[4], lwd = 4, lty = 3)
abline(h = 5, col = cols[5], lwd = 4, lty = 1)
abline(h = 6, col = cols[6], lwd = 4, lty = 1)
abline(h = 7, col = cols[7], lwd = 4, lty = 1)