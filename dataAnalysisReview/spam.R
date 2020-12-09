
# Main simulation function for the sparse additive modeling (SpAM) of
# Ravikumar.

simulation.spam <- function(dat_train, dat_test, dat_val,
                            nbasis = 3, lambda.min.ratio = 1e-3, ...) {
  require(SAM)
  p <- ncol(dat_train) - 1
  n <- nrow(dat_train)

  # Fit the model on training set.
  mod <- samLL(X = dat_train[,-1],
               y = as.numeric(as.character(dat_train$y)), p = nbasis,
               nlambda = 50,
               lambda.min.ratio = lambda.min.ratio, ...)

  spam_res <- eval.obj(mod, dat_train, dat_test, dat_val,
                        pred.fun = predict.spamLL,
                       name = paste0("SpAM", nbasis))
  sparsity_auc <- mean(mod$func_norm[,spam_res[3]] == 0)
  sparsity_mse <- mean(mod$func_norm[,spam_res[4]] == 0)

  return(data.frame("val_auc" = spam_res["auc"],
                    "val_mse" = spam_res["mse"],
                    "sparse_auc" = sparsity_auc,
                    "sparse_mse" = sparsity_mse))

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
  out$values = cbind(Zt, rep(1, nt)) %*% rbind(object$w, object$intercept)
  rm(Zt, newdata)
  return(out)
}

predict.spamLL = function (object, newdata, thol = 0.5) {
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
    Zt[, tmp] = ns(newdata[, j], df = object$p,
                   knots = object$nkots[,j],
                   Boundary.knots = object$Boundary.knots[, j])
  }
  out$probs = exp(cbind(Zt, rep(1, nt)) %*% object$w)
  out$probs = out$probs/(1 + out$probs)
  out$labels = sign(out$probs > thol)
  rm(Zt, newdata)
  return(out$probs)

}
