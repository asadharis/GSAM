
library(glmgen)

###############################################################
#
# This file consists of a proximal gradient descent algorithm
# implementation for the trend-filtering problem.
#
###############################################################

# A function for the proximal operator which solves:
#
# minimize (1/2n) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda1)*||f|| +(lambda2)*TV(f)
#
# where TV(f) is the trend filtering penalty.
# Args:
#   y.ord: The response vector, ordered according to the ordering of x.
#   x.ord: The sorted covariate vector.
#   k: The order of the trend-filtering. We consider 0,1,2
#   lambda1,lambda2: The two tuning parameters.
#   ...: Other parameters that can be passed onto the function glmgen::trendfilter
#
# Returns:
#   f_hat: An n-vector, the solution of the optimization problem.
solve.prox.tf <- function(y.ord, x.ord, k = 0, lambda1, lambda2, ...) {

  # We first solve the problem
  #      f^hat_lambda2 <- argmin (1/2n) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda2)*TV(f)
  # followed by a soft thresholding step.

  require(glmgen)
  n <- length(y.ord)

  # We use n*lambda2 since the default function solves
  #     argmin (1/2) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda2)*TV(f)
  # whereas we would like to have (1/2n) instead of (1/2).

  f_hat <- trendfilter(x = x.ord, y = y.ord,  k = k, family = "gaussian",
                       lambda = n*lambda2, thinning = TRUE, ...)#$beta[, 1]

  # The ... allow for other parameters to be passed, we normally would like to
  # pass the control object, e.g. the control below increases tolerance and the
  # maximum number of iterations.
  #
  # control = trendfilter.control.list(obj_tol = 1e-12, max_iter = 600)

  # The following procedure is used becaused of the thinning feature
  # of trendfiltering solver. While this is very useful for making the solver
  # more efficient this can lead to resulting vectors being smaller than n.
  # The if{} block below aims to rectify this issue.

  if(length(f_hat$beta) != n){
    f_hat <- approx(x = f_hat$x, y = f_hat$beta[,1], xout = x.ord,rule = 2)$y
    # # Find how many values are missing
    # n.res <- n - length(f_hat)
    # # Find where in x.ord we have values too close to each other. And find the smallest 'n.res' values
    # ind.s <- order(diff(x.ord))[1:n.res]
    # f_hat <- R.utils::insert(f_hat, ind.s, values = f_hat[ind.s])
  } else {
    f_hat <- f_hat$beta
  }

  # We return the desired value after soft-thresolding.
  max(1 - lambda1/sqrt(mean(f_hat^2)), 0)*f_hat
}

