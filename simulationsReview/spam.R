library(SAM)

SimSPAM <- function(dat, p = 3,...) {

  fit <- samQL(dat$x, dat$y, p = p, ...)
  #fit <- samQL(dat$x, dat$y, lambda.min.ratio = 1e-4, nlambda = 50, p  = 3)

  # Evaluate the MSE_true
  yhat <- predict.spam(fit, newdata = dat$x)$values
  mse.true <- colMeans((yhat - dat$f0)^2)

  # plot(fit$lambda,mse.true, log = "yx")

  # Find index of best min MSE_test
  yhat.test <- predict.spam(fit, dat$x.test)$values
  mse.test <- colMeans((yhat.test - dat$y.test)^2)

  # Find the index which minimizes the MSE
  ind <- which.min(mse.test)

  # Find  MSE_TRUE_BEST
  mse.true.best <- mse.true[ind]

  # Find the number of active features in best set.
  active.set <- apply(fit$w, 2, function(vec){
    sum(colSums(abs(matrix(vec, nrow = fit$p))) != 0)*1
  })

  # Find the fitted functions of the best lambda value
  # i.e. the lambda value which minizes the test Error.
  fhat.best <- spam.est.func(fit.spam = fit, dat$x, ind)
  xout <- seq(-2.5, 2.5, length = 1000)
  f1 <- approx(dat$x[,1], fhat.best[,1], xout = xout, rule = 2)$y
  f2 <- approx(dat$x[,2], fhat.best[,2], xout = xout, rule = 2)$y
  f3 <- approx(dat$x[,3], fhat.best[,3], xout = xout, rule = 2)$y
  f4 <- approx(dat$x[,4], fhat.best[,4], xout = xout, rule = 2)$y

  fhat.best2 <- cbind(f1,f2,f3,f4)


  return(list("mse.true.best" = mse.true.best,
              "mse.true" = mse.true,
              "lam" = fit$lambda,
              "act.set" = active.set,
              "ind" = ind, "fhat" = fhat.best2))
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


