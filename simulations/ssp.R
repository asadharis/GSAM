library(uniSolve)

SimSPLINE <- function(dat, sq.norm = FALSE, ...) {
  x <- dat$x
  y <- dat$y

  fit <- sobolev.norm(dat$y, dat$x, norm.sq = sq.norm, ...)
  # fit <- sobolev.norm(dat$y, dat$x, norm.sq = sq.norm, lambda.max = 1,
  #                     lambda.min.ratio = 1e-3)

  fit.vals <- apply(fit$f_hat, c(1,3),sum)
  fit.vals <- fit.vals + fit$y.mean

  # Evaluate the MSE_true
  mse.true <- colMeans((fit.vals - dat$f0)^2)

  # Find index of best min MSE_test
  yhat.test <- predict.sobolev(fit, dat$x.test)
  yhat.test <- apply(yhat.test, c(1,3),sum)
  yhat.test <- yhat.test + fit$y.mean

  mse.test <- colMeans((yhat.test - dat$y.test)^2)
  ind <- which.min(mse.test)

  # Find MSE_val and MSE_TRUE_BEST
  yhat.val <- predict.sobolev(fit, dat$x.val)
  yhat.val <- apply(yhat.val, c(1,3),sum)
  yhat.val <- yhat.val + fit$y.mean

  mse.val <- mean((yhat.val[, ind] - dat$y.val)^2)
  mse.true.best <- mean((dat$f0 - fit.vals[, ind])^2)

  # Find the number of active features.
  active.set <- 1*(apply(abs(fit$f_hat), c(2,3), mean)!=0)
  active.set <- colSums(active.set)


  return(list("mse.true.best" = mse.true.best,
              "mse.val" = mse.val,
              "mse.true" = mse.true,
              "lam" = as.numeric(fit$lam),
              "act.set" = active.set,
              "ind" = ind))
}


