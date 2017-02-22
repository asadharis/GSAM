set.seed(1)
n <- 100

x <- rnorm(n)
y <- rnorm(n) + x
theta <- rep(0,length(y))

lambda <- 10

ord <- order(x)
x.ord <- x[ord]
y.ord <- y[ord]

library(uniSolve)
junk <- tf_dp(y,y.ord,lambda,theta)


junk2 <- spline_s(y,y.ord,x.ord, theta,lambda)



# set.seed(1)
# n <- 10
# y <- rnorm(n)
# N <- matrix(rnorm(n*n), ncol = n)
# sigma <- matrix(rnorm(n*n), ncol = n)
# sigma <- sigma %*%t(sigma)
#
#
# myf <- function(lam, y,N,sigma) {
#
#   temp <- t(N)%*% y
#   interm <- crossprod(N) + lam*sigma
#   beta <- solve(interm, temp)
#   lam*as.numeric(crossprod(beta, sigma) %*% beta)
#
# }
#
# xs <- seq(0.1,10,length = 100)
# as <- sapply(xs,myf,y=  y, N = N, sigma = sigma)
# plot(xs,as)
