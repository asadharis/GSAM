# We now solve the optimization problem.
# minimize (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda1)*||f|| +(lambda2)*sqrt(J(f))
#
# To do this we first solve the prox problem
#   f^hat_lambda2 <- argmin (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda2)*sqrt(J(f))
#
# Which requires us to find the root of the equation
# x*sqrt(J)(f^tilde_x) - lambda2/2
# where
# f^tilde_x <- argmin (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (x)*J(f)
#
# Args:
#   y.ord: response vector, n-vector ordered.
#   x.ord: covariate, x-vector also ordered.
#   lambda1, lambda2: Two tuning parameters.
#
# Returns:
#   An n-vector of solution.
#
solve.prox.spline <- function(y.ord, x.ord, lambda1, lambda2) {
  require(stats)
   n <- length(y.ord)
  tempf <- function(lam_x) {
    cpp_temp_func(lam_x, y.ord, x.ord,
                  n, n_grid = 1000, lambda2)
  }
  if(tempf(lambda2*1e+3) < 0) {
    b1 <- cov(x.ord,y.ord)/var(x.ord)
    b0 <- mean(y.ord) - (b1*mean(x.ord))
    fhat <- b0 + (b1*x.ord)
  } else {
    lam <- uniroot(tempf, c(0,lambda2*1e+3),
                   tol = min(lambda2^2,.Machine$double.eps^0.25))$root
    f_hat <- cpp_spline(y.ord, x.ord, lam, n, 1000)
    fhat <- f_hat$sy
  }

  # Return final value.
  max((1 - lambda1/sqrt(mean(fhat^2))), 0)*fhat
}
