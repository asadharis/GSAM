
library(glmnet)

# First a cross-validation for HierBasis
cv.lasso <- function(x.train, y.train, folds, lam.seq, ...){

  require(glmnet)
  n <- length(y.train)
  num.folds <- length(unique(folds))
  nlam <- length(lam.seq)

  pred.errors <- matrix(NA, ncol = nlam, nrow = num.folds)

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

    mod <- glmnet(temp.x, temp.y, family = "gaussian",
                  lambda = lam.seq , ...)

    temp.new.x <- scale(x.train, center = xbar, scale = x.sd)

    preds <- predict(mod, newx = temp.new.x, type = "response")[-temp.ind, ]
    pred.errors[i, ] <- apply((preds - temp.new.y)^2, 2, mean)
  }

  mu.err <- apply(pred.errors, 2, mean)

  se.err <- apply(pred.errors, 2, sd)/sqrt(num.folds)

  return(list("mu" = mu.err, "se" = se.err))
}


simulation.lasso <- function(x.train, y.train, x.test, y.test,
                             folds, lambda.min.ratio, ...) {
  require(glmnet)
  p <- ncol(x.train)
  n <- nrow(x.train)

  cat("Before full model", "\n")
  full.mod <- glmnet(x.train, y.train, family = "gaussian",
                     nlambda = 50, lambda.min.ratio = lambda.min.ratio)
  cat("Full model done", "\n")

  lams <- full.mod$lambda
  cat("Full model done", "\n")
  cv <- cv.lasso(x.train, y.train, folds = folds,
                  lam.seq = lams)

  cat("CV done", "\n")

  # Obtain the minimum CV ONLY!
  ind.min <- which.min(cv$mu)[1]
    ind.1se <- which(cv$mu[ind.min]  - cv$se[ind.min] <= cv$mu &
                       cv$mu[ind.min]  + cv$se[ind.min] >= cv$mu)[1]
  preds <- predict(full.mod, newx = x.test, type = "response")
  ans <- apply((preds[, c(ind.min, ind.1se)] - y.test)^2, 2, mean)
  names(ans) <- c("min", "onese")

  betas <- full.mod$beta[, c(ind.min, ind.1se)]
  act.set <- apply(betas, 2, function(x){
    which(x !=0)
  })
  if(class(act.set) == "matrix") {
    act.set <- list(act.set[,1], act.set[,2])
  }
  if(is.numeric(act.set)) {
    act.set <- list(act.set[1], act.set[2])
  }
  names(act.set) <- names(ans)

  lam.val <- full.mod$lam[c(ind.min, ind.1se)]
  names(lam.val) <- names(ans)

  min.plot <- list("xmat" = x.train[,act.set$min])
  yhat <- t(apply(min.plot$xmat, 1, "*", full.mod$beta[act.set$min, ind.min]))
  min.plot$yhat <- scale(yhat, scale = FALSE)

  onese.plot <- list("xmat" = x.train[,act.set$onese])
  yhat <- t(apply(onese.plot$xmat, 1, "*", full.mod$beta[act.set$onese, ind.1se]))
  onese.plot$yhat <- scale(yhat, scale = FALSE)
  return(list("err" = ans, "sparse" = act.set,
              "lam.val" = lam.val,
              "min.plot" = min.plot,
              "onese.plot" = onese.plot,
              "ind" = c(ind.min, ind.1se)))
}

