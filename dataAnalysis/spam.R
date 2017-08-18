# This file contains helper functions for the method
# SPAM: Sparse Additive Model

library(SAM)

# A cross-validation for SPAM
cv.spam <- function(x.train, y.train, folds, lam.seq, ...){

  #require(SAM)
  n <- length(y.train)
  num.folds <- length(unique(folds))

  pred.errors <- matrix(NA, ncol = length(lam.seq), nrow = num.folds)
  for(i in 1:num.folds) {
    cat("CV fold number: ", i, "\n")
    # Obtain the index of things we will train on this round.
    temp.ind <- which(folds != i)

    # Obtain the response for the covariates and the
    temp.y <- y.train[temp.ind]
    temp.new.y <- y.train[-temp.ind]

    # Scale the hold-out wrt to the training folds.
    temp.x <- scale(x.train[temp.ind, ])
    xbar <- attributes(temp.x)$'scaled:center'
    x.sd <- attributes(temp.x)$'scaled:scale'

    #mod <- samLL(temp.x, temp.y, lambda = lam.seq, p=nbasis)
    mod <- samLL(temp.x, temp.y, lambda = lam.seq, ...)

    temp.new.x <- scale(x.train[-temp.ind,], center = xbar, scale = x.sd)

    preds <- predict.spam(mod, newdata = temp.new.x)$probs
    pred.errors[i, ] <- apply((round(preds) - temp.new.y)^2, 2, mean)
  }

  mu.err <- apply(pred.errors, 2, mean)

  se.err <- apply(pred.errors, 2, sd)/sqrt(num.folds)

  return(list("mu" = mu.err,
              "se" = se.err))

}

simulation.spam <- function(x.train, y.train, x.test, y.test,
                            folds, nbasis = 3, lambda.min.ratio = 1e-4, ...) {
  #require(SAM)
  p <- ncol(x.train)
  J <- nbasis

  cat("Before full model", "\n")
  full.mod <- samLL(x.train, y.train,
                    nlambda = 50, lambda.min.ratio = lambda.min.ratio,
                    p = nbasis, ...)

  cat("Full model done", "\n")
  cv <- cv.spam(x.train, y.train, folds = folds,
                    lam.seq = full.mod$lambda, p = nbasis, ...)
  cat("CV done", "\n")

  # Obtain the minimum CV and one SE lambda index
  ind.min <- which.min(cv$mu)[1]
  ind.1se <- which(cv$mu[ind.min]  - cv$se[ind.min] <= cv$mu &
          cv$mu[ind.min]  + cv$se[ind.min] >= cv$mu)[1]
  preds <- predict.spam(full.mod, newdata = x.test)$probs
  ans <- apply( (round(preds[, c(ind.min, ind.1se)]) - y.test)^2, 2, mean)
  names(ans) <- c("min", "onese")

  # Then we obtain the active set indices.
  betas <- full.mod$w[-(p*J + 1), c(ind.min, ind.1se)]
  act.set <- apply(betas, 2, function(x) {
    temp <- matrix(x, ncol = p, nrow = J)
    ind <- which(colSums(abs(temp))!=0)

    if(length(ind) == 0) {
      return(NaN)
    } else {
      return(ind)
    }
  })
  if(class(act.set) == "matrix") {
    act.set <- list(act.set[,1], act.set[,2])
  }
  if(is.numeric(act.set)) {
    act.set <- list(act.set[1], act.set[2])
  }

  # Then we obtain the lambda values which gave us these results.
  names(act.set) <- names(ans)
  lam.val <- full.mod$lam[c(ind.min, ind.1se)]
  names(lam.val) <- names(ans)

  # Finally we gather information for plotting.
  # For minCV first
  min.yhat <- spam.est.func(full.mod, x.train, index.lam = ind.min)
  min.yhat <- min.yhat[, act.set$min]
  min.xmat <- x.train[, act.set$min]

  # For oneSE second
  onese.yhat <- spam.est.func(full.mod, x.train, index.lam = ind.1se)
  onese.yhat <- onese.yhat[, act.set$onese]
  onese.xmat <- x.train[, act.set$onese]


  return(list("err" = ans, "sparse" = act.set,
              "lam.val" = lam.val,
              "min.plot" = list("yhat" = min.yhat, "xmat" = min.xmat),
              "onese.plot" = list("yhat" = onese.yhat, "xmat" = onese.xmat)
              ))
}



####################################
# Other helper functions for SPAM ##
####################################


# This function obtains the function plots for a
# given lambda index.
# It returns an n by d matrix where d is the number of
# covariates.
spam.est.func = function(fit.spam, X, index.lam) {

  index <- index.lam
  newdata <- X


  #fit.spam$w = fit.spam$w[-nrow(fit.spam$w),]
  # d is number of covariates, fit.spam$p is number of basis functions
  d = ncol(X); n = nrow(X)
  nt = nrow(X)
  X.min.rep = matrix(rep(fit.spam$X.min,nt),nrow=nt,byrow=T)
  X.ran.rep = matrix(rep(fit.spam$X.ran,nt),nrow=nt,byrow=T)

  newdata = (newdata - X.min.rep)/X.ran.rep
  newdata = pmax(newdata, 0)
  newdata = pmin(newdata, 1)
  m = fit.spam$p * d
  ans <- matrix(0, nrow = n, ncol = d)
  for (j in 1:d) {
    tmp = (j - 1) * fit.spam$p + c(1:fit.spam$p)
    tx = ns(newdata[, j], df = fit.spam$p, knots = fit.spam$nkots[,j],
                   Boundary.knots = fit.spam$Boundary.knots[, j])
    ans[,j] <- tx %*% fit.spam$w[tmp, index]
  }
  return(scale(ans, scale = FALSE))
}

predict.spam = function (object, newdata, thol = 0.5) {

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
    tmp = (j - 1) * object$p + c(1:object$p)
    Zt[, tmp] = ns(newdata[, j], df = object$p, knots = object$nkots[,j],
                   Boundary.knots = object$Boundary.knots[, j])
  }
  out$probs = exp(cbind(Zt, rep(1, nt)) %*% object$w)
  out$probs = out$prob/(1 + out$prob)
  out$labels = sign(out$values > thol)
  rm(Zt, newdata)
  return(out)
}


