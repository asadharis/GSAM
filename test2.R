library(GSAM)
set.seed(1)
n <- 300
p <- 8300
x <- matrix(rnorm(n*p), ncol = p, nrow = n)

y <- sin(2*x[,1]) + x[,2]^2 + rnorm(n, sd = 0.1)
library(glmgen)
library(parallel)
library(doParallel)
system.time(mod.tf.paraF <- fit.additive(y, x, max.iter = 100, tol = 1e-5,
                                         family = "gaussian",lambda.max = 4,
                                         lambda.min.ratio = 1e-2,
                                         coord.desc = FALSE,
                                         method = "tf", parallel = FALSE))

system.time(mod.tf.paraT <- fit.additive(y, x, max.iter = 100, tol = 1e-5,
                                         family = "gaussian",lambda.max = 4, lambda.min.ratio = 1e-2,
                                         coord.desc = FALSE, method = "tf",
                                         parallel = TRUE))



myf <- function(t,y) {
  exp(y*t)/((1 + exp(y*t))^2)
}

xs <- seq(-10,10, length = 1e+3)
plot(xs, myf(xs, y=1))
