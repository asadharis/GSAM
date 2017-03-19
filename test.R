set.seed(1)
n <- 100

x <- runif(n, min = -2.5, max = 2.5)
y <- sin(3*x) + rnorm(n, sd = 0.1)
theta <- rep(0,length(y))
par(mfrow = c(1,1))
plot(x,y)
# lambda <- 0.01

ord <- order(x)
x.ord <- x[ord]
y.ord <- y[ord]

obj <- spline_s(y,y.ord, x.ord, theta, 1e-5)

plot(x,y)
lines(obj$x,obj$sy, type = "l", col = "red")

obj2 <- tf_dp(y,y.ord,1,rep(0,n))
lines(x.ord, obj2$beta)

obj3 <- tf_test(y,y.ord,1)
lines(x.ord, obj3$beta, col = "blue")

library(uniSolve)
#junk <- tf_dp(y,y.ord,lambda,theta)

set.seed(1)
x  <- runif(200, -2.5,2.5)
y  <- sin(2*x) + rnorm(200, sd = 0.3)
x <- cbind(x, runif(200, -2.5,2.5),runif(200, -2.5,2.5),runif(200, -2.5,2.5))
y.mean <- mean(y)
x.mean <- apply(x,2,mean)
ord <- apply(x,2,order)
y <- y-mean(y)
x <- scale(x, scale = FALSE)


max.iter <- 100
tol <- 1e-4
initpars <- NULL; lambda1 <- 0.03;
lambda2 <- lambda1^2

ans <- ssp_one(y,y.mean,x,x.mean, ord, lambda1,lambda2)

