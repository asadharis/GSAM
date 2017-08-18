# Main function for a univariate spline solver
spline_s <- function(y, y.ord, x.ord, theta, lambda){
  # According to the code the default function solves the optimization
  # problem.
  # minimize  (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda)*J(f)
  #
  n <- length(x.ord)
  n.grid <- 500
  .C("callSS_Fortran", y = as.double(y.ord),
     x = as.double(x.ord), sy = as.double(theta),
     lambda = as.double(lambda), n_point = as.integer(n),
     yg = as.double(rep(0,n.grid)), xg = as.double(rep(0,n.grid)), ngrid = as.integer(n.grid),
     iderv = as.integer(2))
}

get.pen <- function(obj) {
  n <- length(obj$xg)
  sum(head(obj$yg, -1)^2*(diff(range(obj$xg))/(n+1)))
}


# Solves the proximal problem for the sobolev norm
solve.prox <- function(y, y.ord, x.ord, theta, lambda1, lambda2) {
  # We now want to solve the optimization problem.
  # minimize (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda1)*||f|| +(lambda2)*sqrt(J(f))
  #
  # To do this we first solve the prox problem
  #   f^hat_lambda2 <- argmin (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda2)*sqrt(J(f))
  #
  # Which requires us to find the root of thf_hate equation
  # x*sqrt(J)(f^tilde_x) - lambda2/2
  # where
  # f^tilde_x <- argmin (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (x)*J(f)

  require(stats)

  n <- length(y)
  tempf <- function(lam_x) {
    obj <- spline_s(y, y.ord, x.ord, theta = rep(0, n), lam_x)
    lam_x*sqrt(get.pen(obj)) - lambda2/2
  }
  if(tempf(lambda2*1e+2) < 0) {
    #obj <- spline_s(y, y.ord, x.ord, theta = rep(0, n), lambda2)
    b1 <- cov(x.ord,y.ord)/var(x.ord)
    b0 <- mean(y.ord) - (b1*mean(x.ord))
    f_hat <- b0 + (b1*x.ord)
  } else {
    lam <- uniroot(tempf, c(lambda2*1e-10,lambda2*1e+2))$root
    f_hat <- spline_s(y, y.ord, x.ord, theta = rep(0, n), lam)$sy
  }

  # Finally we return the desired value
  max((1 - lambda1/mean(f_hat^2)), 0)*f_hat

}


ssp_one <- function(y, y.mean, x, x.mean, ord, lambda1, lambda2, max.iter = 100, tol = 1e-4,
                    initpars = NULL) {
  # THis function solves the optimization problem:
  #
  # argmin (1/2n) ||y - sum_{j=1}^{p}f_j||_n^2 +
  #                       \sum_{j=1}^p (lambda1)||f_j||_n +
  #                                    (lambda2)*sqrt(J(f_j)).
  #
  #
  #
  # y: A response vector assumed to be centered
  # y.mean: Mean of uncentered response vector
  # x: A n*p matrix of covariates assumed to be column centered
  # x.mean: Means of uncentered design x.
  # ord: Matrix of dim n*p giving the orders/ranks of each coavairte.
  # lambda1, lambda2: Scalar tuning parameters
  # max.iter: maximum number of iterations for block coordinate descent
  # tol: Tolerance for algorithm
  # intpars: Initial parameters, taken as 0 if none provided

  n <- length(y)
  p <- ncol(x)

  if(is.null(initpars)) {
    initpars <- matrix(0, ncol = p, nrow = n)
  }

  old.pars <- initpars

  # Begin loop for block coordinate descent
  for(i in 1:max.iter) {
    for(j in 1:p) {
      res <- y - apply(initpars[, -j],1,sum)
      initpars[ord[,j], j] <- solve.prox(res, res[ord[,j]], x[ord[,j],j], rep(0, n), lambda1, lambda2)
    }
    if(mean((initpars - old.pars)^2) < tol ) {
      return(initpars)
    } else {
      old.pars <- initpars
    }
  }
  return(initpars)
}



sobolev.norm <- function(y, x, family = "binomial",max.iter = 100, tol = 1e-4,
                         initpars = NULL, initintercept = 0,
                         lambda.max = 1,
                         lambda.min.ratio = 1e-3,
                         nlam = 50, step = 1, alpha = 0.5,
                         gamma.par = 0.5) {
  # This function solves the sobolev norm prooblem
  # In this function we will solve the optmization problem for a
  # sequence of lambda values using warm starts.
  # We solve the above lambda values for the lambda1 = lambda and
  # lambda2 = lambda^2.
  # set.seed(1)
  # n <- 200
  # p <- 4
  # y = rbinom(n, size = 1, prob = 0.5); x = matrix(rnorm(n*p), ncol = p);
  # max.iter = 1e+4; tol = 1e-4;
  # initpars = NULL;
  # initintercept = NULL;
  # lambda.max = 1;
  # lambda.min.ratio = 1e-6;
  # nlam = 50
  # step = 1
  # alpha = 0.5
  # family = "binomial"

  if(family == "binomial"){
    y[y==0] <- -1
  }

  n <- length(y)
  p <- ncol(x)

  if(family == "gaussian") {
    x.mean <- apply(x, 2, mean)
    y.mean <- mean(y)


    x.cen <- scale(x, scale = FALSE)
    y.cen <- y - y.mean
    ord <- apply(x.cen, 2, order)
  } else {
    ord <- apply(x, 2, order) - 1
    ranks <- apply(x, 2, rank) - 1
    x_ord <- apply(x, 2, sort)
  }

  lam.seq <- seq(log10(lambda.max), log10(lambda.max*lambda.min.ratio),
                 length = nlam)
  lam.seq <- 10^lam.seq

  if(is.null(initpars)) {
    initpars <- matrix(0, ncol = p, nrow = n)
  }

  if(is.null(initintercept) & family == "binomial") {
    initintercept = 0;
  }

  ans <- array(0, dim = c(n,p,nlam))
  ans.inters <- numeric(nlam)

  for(i in 1:nlam) {
    if(family == "gaussian") {
      ans[,, i] <- ssp_one(y.cen, y.mean, x.cen, x.mean,
                           ord, lam.seq[i], lam.seq[i],
                           max.iter = max.iter,
                           tol = tol, initpars = initpars)
      initpars <- ans[,, i]
    } else if(family == "binomial") {
      #cat("Num lambda: ", i, "\n")
      temp <- cpp_spp_one(y, x_ord, ord, ranks,
                          lambda1 = lam.seq[i], lambda2 = lam.seq[i]^2,
                          initpars, initintercept, n, p,
                          max_iter = max.iter,  tol = tol,
                          step_size = step, alpha = alpha)
      ans[,,i] <- temp$fhat
      ans.inters[i] <- temp$intercept
      initpars <- temp$fhat
      initintercept <- temp$intercept

    }

  }

  if(family == "gaussian") {
    obj <- list("f_hat" = ans,
                "x.cen" = x.cen,
                "y.cen" = y.cen,
                "x.mean" = x.mean,
                "y.mean" = y.mean,
                "ord" = ord,
                "lam" = lam.seq,
                "family" = family)
  } else if(family == "binomial") {
    obj <- list("f_hat" = ans,
                "intercept" = ans.inters,
                "x" = x,
                "ord" = ord,
                "lam" = lam.seq,
                "family" = family)
  }

  class(obj) <- "sobolev"
  return(obj)
}


predict.sobolev <- function(obj, new.data, type = "function") {

  # new.data <- matrix(rnorm(n*p), ncol = p)
  if(obj$family == "gaussian") {
    new.dat.cen <- scale(new.data, scale = FALSE, center = obj$x.mean)

    nlam <- length(obj$lam)
    p <- dim(obj$f_hat)[2]

    ans <- array(0, dim = c(dim(new.data), nlam) )
    for(i in 1:nlam) {
      for(j in 1:p) {
        ans[,j,i] <- approx(obj$x.cen[, j], obj$f_hat[,j,i], new.dat.cen[, j],
                            rule = 2)$y
      }
    }
  } else if(obj$family == "binomial") {

    nlam <- length(obj$lam)
    p <- dim(obj$f_hat)[2]

    ans <- array(0, dim = c(dim(new.data), nlam) )
    for(i in 1:nlam) {
      for(j in 1:p) {
        ans[,j,i] <- approx(obj$x[, j], obj$f_hat[,j,i], new.data[, j],
                            rule = 2)$y
      }
    }
  }
  if(type == "function") {
    return(ans)
  } else {
    if(obj$family == "gaussian"){
      return(apply(ans, c(1,3),sum) + obj$y.mean)
    } else {
      temp <- apply(ans, c(1,3),sum)
      temp <- t(apply(temp, 1, "+", obj$intercept))

      return(1/(1 + exp(-temp)))

    }
  }

}

