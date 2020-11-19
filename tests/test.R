set.seed(1)
n <- 500
x.ord <- seq(-pi, pi, length = n)
y.ord <- sin(2*x.ord) + rnorm(n, sd = 0.2)

plot(x.ord, y.ord)



library(GSAM)
lam <- 0.05
lambda2 <- lam^2
tempf <- function(lam_x) {
  GSAM:::cpp_temp_func(lam_x, y.ord, x.ord,
                       n, n_grid = 1000, lambda2)
}

library(microbenchmark)

r1 <- function() {
  uniroot(tempf, c(0,lambda2*1e+3),
          tol = min(lambda2^2,1e-10))$root
}

c1 <- function() {
  GSAM:::cpp_uniroot(0, lambda2*1e+3,y.ord,x.ord,
                     n, n_grid=1000, lambda2,
                     tol = min(lambda2^2, 1e-10))
}
c2 <- function() {
  GSAM:::cpp_uniroot2(0, lambda2*1e+3,y.ord,x.ord,
                      n, n_grid=1000, lambda2,
                      tol = min(lambda2^2, 1e-10),
                      max_iter = 3000)
}
c3 <- function() {
  GSAM:::cpp_uniroot3(0, lambda2*1e+3, min(lambda2^2, 1e-10), y.ord, x.ord,n,
                      n_grid = 1000, lambda2 = lambda2)
}

microbenchmark(r1(),c1(),c2(), c3())


microbenchmark( a <- solve.prox.spline(y.ord, x.ord, lam, lam^2),
  b <- cpp_solve_prox(y.ord, x.ord, lam, lam^2, n, 1000, 1)
)



lines(x.ord, a, col  = "red", lwd = 2)
lines(x.ord, b, col  = "blue", lwd = 2)

#################################################################
#################################################################
## TESTING time difference between splines and TF
#################################################################
#################################################################

library(GSAM)
set.seed(1)
n <- 500
x.ord <- seq(-pi, pi, length = n)
y.ord <- sin(2*x.ord) + rnorm(n, sd = 0.2)

plot(x.ord, y.ord)
lambda1 <-  0.001
lambda2 <- lambda1^2

times <- microbenchmark(
  "sp" = sp <- cpp_solve_prox(y.ord, x.ord, lambda1, lambda2, n, 1e+3, 0.5),
  "tf0" = tf0 <- solve.prox.tf(y.ord, x.ord, k=0, lambda1, lambda2),
  "tf1" = tf1 <- solve.prox.tf(y.ord, x.ord, k=1, lambda1, lambda2),
  "tf2" = tf2 <- solve.prox.tf(y.ord, x.ord, k=2, lambda1, lambda2)
  )

boxplot(times)
# It looks like adding a bisection to splines makes it slower than tf0
# but still faster than tf1 and tf2. This is encouraging as tf0 is fused lasso
# with known algorithms that are super fast.

#################################################################
#################################################################
#################################################################

set.seed(1)
n <- 300
p <- 83
x <- matrix(rnorm(n*p), ncol = p, nrow = n)

y <- sin(2*x[,1]) + x[,2]^2 + rnorm(n, sd = 0.1)
f_hat <- cbind(sin(2*x[,1]), x[,2]^2, x[,-(1:2)])
intercept <- mean(y)
step_size <- 0.001

# This way is faster than the current implementation in package.
ord_mat <- apply(x,2,order)
x_mat_ord <- sapply(1:ncol(ord_mat), function(i){
  x[ord_mat[,i], i]
})

lambda1 <- 0.01
lambda2 <- lambda1^2
y <- y
family = "gaussian"
method = "tf"
k <- 2


ncores <- 8
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

myf1 <- function(){
  GSAM:::GetZ(f_hat, intercept, step_size, x_mat_ord, ord_mat,
            k=k, lambda1, lambda2, y, family, method, parallel = FALSE)
}

myf2 <- function(){
  GSAM:::GetZ(f_hat, intercept, step_size, x_mat_ord, ord_mat,
              k=k, lambda1, lambda2, y, family, method, parallel = TRUE,
              ncores = 8)
}

library(rbenchmark)
benchmark(myf1(), myf2(), replications = 50)

a <- myf1()[[2]]
b <- myf2()[[2]]
parallel::stopCluster(cl)

library(microbenchmark)
microbenchmark(apply(f_hat,1,sum), rowSums(f_hat))

#################################################################
#################################################################
#################################################################
# Final testing of parallelization.
# TESTING package version 0.4.1

library(GSAM)
set.seed(1)
n <- 300
p <- 43
x <- matrix(rnorm(n*p), ncol = p, nrow = n)

y <- sin(2*x[,1]) + x[,2]^2 + rnorm(n, sd = 0.1)
system.time(mod.tf.paraF <- fit.additive(y, x, max.iter = 100, tol = 1e-5,
             family = "gaussian",lambda.max = 4, lambda.min.ratio = 1e-2,
             coord.desc = FALSE, method = "tf"))

system.time(mod.tf.paraT <- fit.additive(y, x, max.iter = 100, tol = 1e-5,
                                         family = "gaussian",lambda.max = 4, lambda.min.ratio = 1e-2,
                                         coord.desc = FALSE, method = "tf",
                                         parallel = TRUE))

#################################################################
#################################################################
#################################################################



