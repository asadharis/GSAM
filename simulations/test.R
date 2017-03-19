
library(glmgen)

set.seed(1)
n = 300
x <- runif(n, -2.5,2.5)
ord <- order(x)
y <- sin(2*x) + rnorm(n, sd = 0.3)

plot(x,y)
glmgenmod <- glmgen::trendfilter(x, y, k = 3, lambda = 0.1)
lines(x[ord],glmgenmod$beta, col = "red")



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
