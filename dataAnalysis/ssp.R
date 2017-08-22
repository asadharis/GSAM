# This file contains all the functions we will use to fit the processed data.

#Rcpp::sourceCpp('helper_functions.cpp')

# First a cross-validation for HierBasis
cv.ssp <- function(x.train, y.train, folds, nlam = 30,
                         max.lambda, lam.min.ratio, gamma.par){
  # We assume that the y.train is coded as (-1,1) for HierBasis.

  require(uniSolve)
  n <- length(y.train)
  num.folds <- length(unique(folds))

  pred.errors <- matrix(NA, ncol = nlam, nrow = num.folds)
  lam.seq <- NULL
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

    mod <- sobolev.norm(temp.y, temp.x, nlam = nlam,
                        lambda.max = max.lambda,
                        lambda.min.ratio = lam.min.ratio,
                        max.iter = 1000, family = "binomial", tol = 1e-8,
                        gamma.par = gamma.par)

    temp.new.x <- scale(x.train[-temp.ind,], center = xbar, scale = x.sd)

    preds <- predict(mod, new.data = temp.new.x, type = "response")
    #preds <- apply(preds, c(1,3),sum) + mod$y.mean
    pred.errors[i, ] <- apply((round(preds) - temp.new.y)^2, 2, mean)

  }

  lam.seq <- as.numeric(mod$lam)

  mu.err <- apply(pred.errors, 2, mean)

  se.err <- apply(pred.errors, 2, sd)/sqrt(num.folds)

  return(list("mu" = mu.err, "se" = se.err, "lam" = lam.seq))
}



simulation.ssp <- function(x.train, y.train, x.test, y.test,
                          folds, max.lambda, lam.min.ratio, gamma.par,...) {
  require(uniSolve)
  p <- ncol(x.train)
  n <- length(y.train)
  cat("Before full model", "\n")

  max.lambda <- 1
  lam.min.ratio <- 0.1
  gamma.par <- NULL
  full.mod <- sobolev.norm(y.train, x.train, lambda.max = max.lambda,
                           lambda.min.ratio = lam.min.ratio, max.iter = 500,
                           tol = 1e-4, gamma.par = gamma.par,step = n,
                           nlam = 10)
  cat("Full model done", "\n")

  apply(abs(full.mod$f_hat),3,function(x){sum(colSums(x)!=0)})

  cv <- cv.ssp(x.train, y.train, folds = folds, max.lambda = full.mod$lam[1],
                     nlam = length(full.mod$lam), lam.min.ratio = lam.min.ratio,
               gamma.par = gamma.par)
  cat("CV done", "\n")
  # Obtain the minimum CV and one SE lambda
  ind.min <- which.min(cv$mu)[1]
  ind.1se <- which(cv$mu[ind.min]  - cv$se[ind.min] <= cv$mu &
                     cv$mu[ind.min]  + cv$se[ind.min] >= cv$mu)[1]

  preds.fhat <- predict(full.mod, new.data = x.test, type = "response")
  ans <- apply((round(preds.fhat[,c(ind.min, ind.1se)]) - y.test)^2, 2, mean)
  names(ans) <- c("min", "onese")

  betas <- full.mod$f_hat[,,c(ind.min, ind.1se)]
  act.min <- which(colSums(abs(betas[,,1])) != 0)
  act.1se <- which(colSums(abs(betas[,,2])) != 0)

  act.set <- list(act.min, act.1se)
  names(act.set) <- names(ans)


  lam.val <- full.mod$lam[c(ind.min, ind.1se)]
  names(lam.val) <- names(ans)

  min.plot <- list("yhat" = full.mod$f_hat[,act.set$min,ind.min],
                   "xmat" = x.train[, act.set$min])
  onese.plot <- list("yhat" = full.mod$f_hat[,act.set$onese,ind.1se],
                     "xmat" = x.train[, act.set$onese])
  return(list("err" = ans, "sparse" = act.set,
              "lam.val" = lam.val,
              "min.plot" = min.plot,
              "onese.plot" = onese.plot))
}
