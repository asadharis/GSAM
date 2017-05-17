# IN this file we have the main functions for the trend-filtering additive models
# We use the c implementation of the package glmgen.

library(glmgen)

# Solves the proximal problem for the sobolev norm
solve.prox.tf <- function(y.ord, x.ord, k = 0, lambda1, lambda2) {
  # We now want to solve the optimization problem.
  # minimize (1/2n) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda1)*||f|| +(lambda2)*TV(f)
  #
  # To do this we first solve the problem
  #   f^hat_lambda2 <- argmin (1/2n) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda2)*TV(f)
  #
  # The x and y are assumed to be centered here.

  require(glmgen)

  n <- length(y.ord)
  f_hat <- trendfilter(x = x.ord, y = y.ord,  k = k, family = "gaussian",
              lambda = n*lambda2)$beta[, 1]

  ######################################################################
  ################ ISSUE IN package glmgen #############################
  ######################################################################
  # THere seems to be an issue in the package glmgen, where it returns solutions
  # of size < n. On inspection, it appears this happens when the x.ord has
  # repeated measurements or measurements <1e-6 apart. The lines of code below
  # aims to resolve this issue.

  if(length(f_hat) != n){
    # Find where in x.ord, two consecutive values are < 1e-6
    ind.s <- which(diff(x.ord) < 1e-6)
    f_hat <- R.utils::insert(f_hat, ind.s, values = f_hat[ind.s])
  }


  #f_hat2 <- genlasso::trendfilter(y = y.ord, pos = x.ord, ord = 0)

  # We return the desired value
  max((1 - lambda1/mean(f_hat^2)), 0)*f_hat
}

tf_one <- function(y, y.mean, x, x.mean, ord, k = 0, lambda1, lambda2,
                       max.iter = 100, tol = 1e-4, initpars = NULL) {
  # THis function solves the optimization problem:
  #
  # argmin (1/2n) ||y - sum_{j=1}^{p}f_j||_n^2 +
  #                       \sum_{j=1}^p (lambda1)||f_j||_n +
  #                                    (lambda2)*TV(f_j).
  #
  #
  # y: A response vector assumed to be centered
  # y.mean: Mean of uncentered response vector
  # x: A n*p matrix of covariates assumed to be column centered
  # x.mean: Means of uncentered design x.
  # ord: Matrix of dim n*p giving the orders/ranks of each coavairte.
  # k: Order of the trend filtering problem. Possible choices are 0,1,2,3.
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

  #p <- 96
  # Begin loop for block coordinate descent
  for(i in 1:max.iter) {
    for(j in 1:p) {
      #print(j)
      res <- y - apply(initpars[, -j], 1, sum)
      initpars[ord[,j], j] <- solve.prox.tf(res[ord[, j]],
                                            x[ord[, j], j], k = k, lambda1, lambda2)
    }

    if(mean((initpars - old.pars)^2) < tol ) {
      return(initpars)
    } else {
      old.pars <- initpars
    }
  }
  return(initpars)
}


tf.norm <- function(y, x, max.iter = 100, tol = 1e-4,
                         initpars = NULL,
                         lambda.max = 1,
                         lambda.min.ratio = 1e-3,
                         nlam = 50,
                         k = 0) {
  # This function solves the trenfiltering/total variation problem.
  # In this function we will solve the optmization problem for a
  # sequence of lambda values using warm starts.
  # We solve the above lambda values for the lambda1 = lambda and
  # lambda2 = lambda^2.

  # n <- 200
  # p <- 100
  # x = matrix(rnorm(n*p), ncol = p);
  # y = sin(3*x[,1]) + rnorm(n,sd = 0.1);
  # max.iter = 100; tol = 1e-4;
  # initpars = NULL;
  # lambda.max = 1;
  # lambda.min.ratio = 1e-3;
  # nlam = 50
  # k <- 0

  n <- length(y)
  p <- ncol(x)

  x.mean <- apply(x, 2, mean)
  y.mean <- mean(y)

  x.cen <- scale(x, scale = FALSE)
  y.cen <- y - y.mean

  ord <- apply(x.cen, 2, order)

  lam.seq <- seq(log10(lambda.max), log10(lambda.max*lambda.min.ratio),
                 length = nlam)
  lam.seq <- 10^lam.seq

  if(is.null(initpars)) {
    initpars <- matrix(0, ncol = p, nrow = n)
  }

  ans <- array(0, dim = c(n,p,nlam))
  for(i in 1:nlam) {
    #print(i)
    ans[,, i] <- tf_one(y.cen, y.mean, x.cen, x.mean,
                            ord, k= k, lam.seq[i], lam.seq[i]^2,
                            max.iter = max.iter,
                            tol = tol, initpars = initpars)
    initpars <- ans[,, i]

  }

  obj <- list("f_hat" = ans,
              "x.cen" = x.cen,
              "y.cen" = y.cen,
              "x.mean" = x.mean,
              "y.mean" = y.mean,
              "ord" = ord,
              "lam" = lam.seq,
              "k" = k)
  class(obj) <- "tf"
  return(obj)
}


predict.tf <- function(obj, new.data) {
  new.dat.cen <- scale(new.data, scale = FALSE, center = obj$x.mean)

  nlam <- length(obj$lam)
  p <- dim(obj$f_hat)[2]

  ans <- array(0, dim = c(dim(new.data), nlam) )
  for(i in 1:nlam) {
    for(j in 1:p) {
      ans[,j,i] <- approx(new.dat.cen[, j], obj$f_hat[,j,i], new.dat.cen[, j])$y
    }
  }
  return(ans)
}


