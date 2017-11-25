
library(glmgen)

set.seed(1)
n = 300
x <- runif(n, -2.5,2.5)
ord <- order(x)
y <- sin(2*x) + rnorm(n, sd = 0.3)

x <- dat$x[,1]
y <- dat$y
ord <- order(x)
plot(x,y)

fhat <-solve.prox.spline(y2[ord[,1]], dat$x[ord[,1],1],
                         fit$lam[32],fit$lam[32]^2)

lines(dat$x[ord[,1],1],fhat, col = "red")

plot(dat$x[ord[,1],1],fhat)

# set.seed(1)
# x  <- runif(200, -2.5,2.5)
# y  <- sin(2*x) + rnorm(200, sd = 0.3)
# x <- cbind(x, runif(200, -2.5,2.5),runif(200, -2.5,2.5),runif(200, -2.5,2.5))
# y.mean <- mean(y)
# x.mean <- apply(x,2,mean)
# ord <- apply(x,2,order)
# y <- y-mean(y)
# x <- scale(x, scale = FALSE)
# max.iter <- 100
# tol <- 1e-4
# initpars <- NULL; lambda1 <- 0.03;
# lambda2 <- lambda1^2
#
# a <- spline_s(y,y[ord[,1]],x[ord[,1],1], rep(0,200), 0.001)
# plot(x[,1],y)
# lines(a$x, a$sy)
