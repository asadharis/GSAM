library(SAM)

SimSPAM <- function(dat, p = 3,...) {
  x <- dat$x
  y <- dat$y

  fit <- samQL(x, y, p = p, ...)
  #fit <- samQL(x, y, lambda.min.ratio = 1e-4, nlambda = 50, p  = 3)

  # Evaluate the MSE_true
  yhat <- predict.spam(fit, newdata = x)$values
  mse.true <- colMeans((yhat - dat$f0)^2)

  # plot(fit$lambda,mse.true, log = "yx")

  # Find index of best min MSE_test
  yhat.test <- predict.spam(fit, dat$x.test)$values
  mse.test <- colMeans((yhat.test - dat$y.test)^2)
  ind <- which.min(mse.test)

  # plot(fit$lambda, mse.test, log = "xy")

  # Find MSE_val and MSE_TRUE_BEST
  yhat.val <- predict.spam(fit, dat$x.val)$values
  mse.val <- colMeans((yhat.val - dat$y.val)^2)[ind]
  mse.true.best <- colMeans((dat$f0 - yhat)^2)[ind]

  # Find the number of active features in best set.
  active.set <- apply(fit$w, 2, function(vec){
    sum(colSums(abs(matrix(vec, nrow = fit$p))) != 0)*1
  })

  return(list("mse.true.best" = mse.true.best,
              "mse.val" = mse.val,
              "mse.true" = mse.true,
              "lam" = fit$lambda,
              "act.set" = active.set,
              "ind" = ind))
}


predict.spam = function (object, newdata){

  gcinfo(FALSE)
  out = list()
  nt = nrow(newdata)
  d = ncol(newdata)
  X.min.rep = matrix(rep(object$X.min, nt), nrow = nt, byrow = T)
  X.ran.rep = matrix(rep(object$X.ran, nt), nrow = nt, byrow = T)
  newdata = (newdata - X.min.rep)/X.ran.rep
  newdata = pmax(newdata, 0)
  newdata = pmin(newdata, 1)
  m = object$p * d
  Zt = matrix(0, nt, m)
  for (j in 1:d) {
    # We fix a typo here! Original version given in comments.
    # The typo is that object$knots[, j] should be object$nkots[, j]
    ############## OLD VERSION #############################
#     tmp = (j - 1) * object$p + c(1:object$p)
#     Zt[, tmp] = ns(newdata[, j], df = object$p,
#                    knots = object$knots[, j], Boundary.knots = object$Boundary.knots[, j])
    ############## OLD VERSION #############################

    tmp = (j - 1) * object$p + c(1:object$p)
    Zt[, tmp] = ns(newdata[, j], df = object$p,
                   knots = object$nkots[, j], Boundary.knots = object$Boundary.knots[, j])
  }
  out$values = cbind(Zt, rep(1, nt)) %*% rbind(object$w, object$intercept)
  rm(Zt, newdata)
  return(out)
}


# A HELPER FOR PLOTTING SPAM
# This function obtains the function plots for a
# given lambda index.
# It returns an n by d matrix where d is the number of
# covariates.
spam.est.func = function(fit.spam, X, index.lam) {

  index <- index.lam

  #fit.spam$w = fit.spam$w[-nrow(fit.spam$w),]
  # d is number of covariates, fit.spam$p is number of basis functions
  d = ncol(X); n = nrow(X)
  nt = nrow(X)
  X.min.rep = matrix(rep(fit.spam$X.min,nt),nrow=nt,byrow=T)
  X.ran.rep = matrix(rep(fit.spam$X.ran,nt),nrow=nt,byrow=T)
  Xt <-  (X - X.min.rep)/(X.ran.rep)
  m = fit.spam$p * d
  ans <- matrix(0, nrow = n, ncol = d)
  for (j in 1:d) {
    tmp = (j-1) * fit.spam$p + c(1:fit.spam$p)
    tx <- ns(Xt[, j], knots = fit.spam$nkots[, j],
             Boundary.knots = fit.spam$Boundary.knots[, j])
    ans[,j] <- tx %*% fit.spam$w[tmp, index]
  }
  return(scale(ans, scale = FALSE))
}


