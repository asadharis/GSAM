
# Main simulation function for the lasso.

simulation.lasso <- function(dat_train, dat_test, dat_val,
                             lambda.min.ratio = 1e-4, ...) {
  require(glmnet)
  p <- ncol(dat_train) - 1
  n <- nrow(dat_train)

  # Fit the model on training set.
  mod <- glmnet(x = as.matrix(dat_train[,-1]), y = dat_train$y,
                family = "binomial",
                nlambda = 50,
                lambda.min.ratio = lambda.min.ratio, ...)

  lasso_res <- eval.obj(mod, dat_train, dat_test, dat_val,
                        pred.fun = pred.lasso, name = "lasso")
  sparsity_auc <- mean(mod$beta[,lasso_res[3]] == 0)
  sparsity_mse <- mean(mod$beta[,lasso_res[4]] == 0)

  return(data.frame("val_auc" = lasso_res["auc"],
              "val_mse" = lasso_res["mse"],
              "sparse_auc" = sparsity_auc,
              "sparse_mse" = sparsity_mse))
}

pred.lasso <- function(obj, dat) {
  predict(obj, newx = dat, type = "response")
}

