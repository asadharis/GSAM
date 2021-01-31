
source("spamProx.R")
# A file to write a wrapper function for all
# of the prox problems.

# Args:
#   r: The output vector, orderd according to x.
#   x: the input vector of ordered x values, can be null if other
#       input provided, e.g. Q matrix for SpAM.
#   lam1: lambda for sparsity penalty.
#   lam2: lambda for smoothness penalty.
#   method: Name of method.
#   k: order for trend filtering.
#   df: number of basis functions for SpAM.
#   ...: Other arguments.
prox_wrapper <- function(r, x = NULL, lam1 = 1, lam2 = 1,
                         method = "spam", k = 0,
                         df = 3, ...) {
  require(GSAM)
  require(glmgen)
  if(method == "spam") {
    prox_spam(r, lam1 = lam1, ...)
  } else if(method == "tf") {
    GSAM::solve.prox.tf(r, x, k = k, lambda1 = lam1,
                        lambda2 = lam2)
  } else if(method == "ssp") {
    n <- length(r)
    GSAM::cpp_solve_prox(r, x, lambda1 = lam1, lambda2 = lam2,
                         n = n,n_grid = 1e+3, lam_tilde_old = 0.5)
  }

}

# ############ TESTING ########################
# # Generate data
# set.seed(1)
# n <- 500
# df <- 30
# x <- sort(runif(n, min = -1,max = 1))
# r <- 0.5*sin(3*x) + rnorm(n, sd = 0.1)
#
# x_matQ <- Qmat_spam(x, df = df)
#
# lam1 =0.04
# lam2 = lam1^2
#
# fit_spam <- prox_wrapper(r, x = NULL, lam1 = lam1, lam2 = lam2,
#                     method = "spam", x_matQ = x_matQ)
# fit_ssp <- prox_wrapper(r, x = x, lam1 = lam1, lam2 = lam2,
#                         method = "ssp")
# fit_tf0 <- prox_wrapper(r, x = x, lam1 = lam1, lam2 = lam2,
#                         method = "tf", k = 0)
# plot(x,r)
# lines(x, fit_spam, col = "red", lwd = 2)
# lines(x, fit_ssp, col = "blue", lwd = 2)
# lines(x, fit_tf0, col = "green", lwd = 2)
