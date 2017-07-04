set.seed(1)
n <- 100

x <- runif(n, min = -2.5, max = 2.5)
pr <- 1/(1+ exp(-(sin(3*x))))
y <- rbinom(n, size = 1, prob = pr)
y[y==0] <- -1
theta <- rep(0,length(y))
par(mfrow = c(1,1))
plot(x,y)
# lambda <- 0.01


x_mat <- cbind(x, runif(n,-2.5,2.5),runif(n,-2.5,2.5),runif(n,-2.5,2.5),runif(n,-2.5,2.5))
ord_mat <- apply(x_mat, 2, function(vec){
  order(vec) - 1
})

rank_mat <- apply(x_mat, 2, function(vec){
  rank(vec) - 1
})


x_mat_ord <- apply(x_mat, 2, function(vec){
  vec[order(vec)]
})

p <- ncol(x_mat)

ans <- cpp_spp_one(y, x_mat_ord,
            ord_mat, rank_mat, lambda1 = 1.1e-3, lambda2 = (1.1e-3)^2,
            init_fhat = matrix(0, ncol = p, nrow = n), init_intercept = 0,
            n, p,
            max_iter = 100000, tol = 1e-4,
            step_size = 1, alpha = 0.4)


plot(x_mat[,1], ans[,1])
plot(x_mat[,2], ans[,2])


ord <- order(x)
x.ord <- x[ord]
y.ord <- y[ord]


ans <- solve.prox(y,y.ord, x.ord, theta,lambda1 = 0.1, lambda2 = 0.1^2)
ans2 <- cpp_solve_prox(y.ord,x.ord,0.1,0.1^2, n, 500)

obj <- spline_s(y,y.ord, x.ord, theta, 1e-5)

plot(x,y)
lines(obj$x,obj$sy, type = "l", col = "red")


obj2 <- cpp_spline(y.ord, x.ord, lambda = 1e-5, n, n_grid = 500)
lines(obj2$x,obj2$sy, type = "l", col = "blue")

obj3 <- cpp_spline_raw(y.ord, x.ord, lambda = 1e-5, n, n_grid = 500)

lambda2 <- 0.01

tempf <- function(lam_x) {
  obj <- spline_s(y, y.ord, x.ord, theta = rep(0, n), lam_x)
  #get.pen(obj)
  lam_x*sqrt(get.pen(obj)) - lambda2/2
}

tempf2 <- function(lam_x) {
  cpp_temp_func(lam_x,y.ord,x.ord,n,500,lambda2)
}

##########################


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

