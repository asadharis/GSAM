# Main simulation function for GSAM with Trend Filtering penalty.

simulation.tf <- function(dat_train, dat_test, dat_val,
                           k = 0, max.lambda = 1,
                           lambda.min.ratio = 1e-3, ...) {
  require(GSAM)
  p <- ncol(dat_train) - 1
  n <- nrow(dat_train)

  # Fit the model on training set.
  mod <- fit.additive(y=dat_train$y, x=dat_train[,-1],
                      family="binomial", method = "tf",
                      k=k,
                      lambda.max = max.lambda,
                      lambda.min.ratio = lambda.min.ratio,
                      max.iter = 100, ...)

  # Evaluate results
  tf_res <- eval.obj(mod, dat_train, dat_test, dat_val,
                      pred.fun = gsam_pred,
                      name = paste0("TF", k))
  # Calculate sparsity
  sparsity_auc <- mean(colMeans(mod$f_hat[, ,tf_res[3]]^2) == 0)
  sparsity_mse <- mean(colMeans(mod$f_hat[, ,tf_res[4]]^2) == 0)

  return(data.frame("val_auc" = tf_res["auc"],
                    "val_mse" = tf_res["mse"],
                    "sparse_auc" = sparsity_auc,
                    "sparse_mse" = sparsity_mse))

}

