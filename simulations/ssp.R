
SimSPLINE <- function(dat, ...) {
  require(uniSolve)

  # fit <- fit.additive(y=dat$y, x=dat$x,
  #                     family="gaussian", method = "sobolev",
  #                     lambda.max = 3, lambda.min.ratio = 1e-2,
  #                     tol = 1e-5, max.iter = 3000)
  fit <- fit.additive(y=dat$y, x=dat$x,
                      family="gaussian", method = "sobolev",
                      ...)

  fit.vals <- apply(fit$f_hat, c(1,3),sum)
  fit.vals <- fit.vals + fit$intercept

  # Evaluate the MSE_true
  mse.true <- colMeans((fit.vals - dat$f0)^2)

  # Find index of best min MSE_test
  yhat.test <- predict(fit, dat$x.test)
  yhat.test <- apply(yhat.test, c(1,3),sum)
  yhat.test <- yhat.test + fit$intercept

  mse.test <- colMeans((yhat.test - dat$y.test)^2)
  ind <- which.min(mse.test)

  # Find MSE_val and MSE_TRUE_BEST
  mse.true.best <- mse.true[ind]

  # Find the number of active features.
  active.set <- 1*(apply(abs(fit$f_hat), c(2,3), mean)!=0)
  active.set <- colSums(active.set)

  # Find the fitted functions of the best lambda value
  # i.e. the lambda value which minizes the test Error.
  fhat.best <- fit$f_hat[, , ind]

  xout <- seq(-2.5, 2.5, length = 1000)
  f1 <- approx(dat$x[,1], fhat.best[,1], xout = xout, rule = 2)$y
  f2 <- approx(dat$x[,2], fhat.best[,2], xout = xout, rule = 2)$y
  f3 <- approx(dat$x[,3], fhat.best[,3], xout = xout, rule = 2)$y
  f4 <- approx(dat$x[,4], fhat.best[,4], xout = xout, rule = 2)$y

  #plot(xout,f4)

  fhat.best2 <- cbind(f1,f2,f3,f4)

  return(list("mse.true.best" = mse.true.best,
              "mse.true" = mse.true,
              "lam" = as.numeric(fit$lam),
              "act.set" = active.set,
              "ind" = ind, "fhat" = fhat.best2))
}


