# Test file for path lemma
library(uniSolve)

set.seed(1)
n <- 500
x <- sort(runif(500,-2.5,2.5))
y <- sin(2*x) + rnorm(n, sd = 0.2)

lam.seq1 <- 10^seq(log10(0.2), log10(0.001), length = 100)
lam.seq <- 10^seq(log10(10), log10(0.00008), length = 100)


spline.ph <- sapply(lam.seq, function(lam){
  solve.spline(y,x,lam = lam)$phat
})

lam.seq2 <- 2 * lam.seq * spline.ph
#plot(lam.seq2, lam.seq, log = "y")


prox.ph <- sapply(lam.seq2, function(lam){
  solve.prox.spline2(y,x,lambda1 = 0, lambda2 = lam)$phat
})
prox.ph2 <- sapply(lam.seq1, function(lam){
  solve.prox.spline2(y,x,lambda1 = 0, lambda2 = lam)$phat
})

plot(lam.seq1, prox.ph2,
     ylim = range(prox.ph2), col = "blue",
     xlab = "",ylab = "", pch = 1, cex = 1,type = "p" )

lines(lam.seq2, spline.ph, col = "red", lwd = 3)

labs <- approx(lam.seq2, lam.seq, xout = c(0, 0.03, 0.08, 0.1167), rule = 2)$y
labs <- round(labs, 3)
labs[4] <- expression(infinity)
axis(side=3, at = c(0, 0.03, 0.08, 0.1167),
     labels =  labs)
mtext(expression(lambda), side = 1, line = 3, cex = 2, col = "blue")
mtext(expression(gamma), side = 3, line = 1.8, cex = 2, col = "red")
mtext(expression(displaystyle(P[st](f))), side = 2, line = 2, cex = 1.5)


axis(side=2)
legend("topright", c("Smoothing Spline", "Proximal Problem"),
       col = c("blue", "red"), lty = c(NA,1),
       pch = c(1,NA), lwd = 2, cex = 1.3)

plot(1:10,1:10)
