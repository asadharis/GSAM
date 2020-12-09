# File for helper functions
#
# Evaluates Method performance.
#
# Input:
#       obj: the fitted method object
#       dat_train/test/val: Datasets with first column 'y', the response.
#       pred.fun: a predict function where first element is object, second is
#                 matrix of predictor values.
#       name: (optional) string for method name (for reading warnings).
# Output: A named vector with
#       auc: AUC on the validation set (estimated AUC) of method.
#       mse: MSE on the validation set (estimated MSE) of method.
#       ind_*: Index of selected method based on *.
eval.obj <- function(obj, dat_train, dat_test, dat_val,
                     pred.fun = predict.spamLL,
                     name = "lasso") {

  test_phat <- pred.fun(obj, as.matrix(dat_test[,-1]))
  val_phat <- pred.fun(obj, as.matrix(dat_val[,-1]))

  nlam <- ncol(test_phat)
  auc_test <- mse_test <- numeric(nlam)
  for(i in 1:nlam){
    #print(i)
    auc_test[i] <- pROC::roc(response = dat_test$y,
                           predictor = test_phat[,i], quiet = TRUE)$auc
    bias <- test_phat[,i] - as.numeric(as.character(dat_test$y))
    mse_test[i] <- mean(bias^2)
  }

  ind_auc <- which.max(auc_test)
  ind_mse <- which.min(mse_test)
  if(ind_auc == 1 | ind_auc == nlam) {
    warning(paste0(name, ": test AUC not concave."))
  }
  if(ind_mse == 1 | ind_mse == nlam) {
    warning(paste0(name, ": test MSE not convex."))
  }

  val_auc <- pROC::roc(response = dat_val$y,
                       predictor = val_phat[,ind_auc],
                       quiet = TRUE)$auc
  val_bias <- val_phat[,ind_mse]- as.numeric(as.character(dat_val$y))
  val_mse <- mean(val_bias^2)

  return(c("auc" = val_auc, "mse" = val_mse,
           "ind_auc" = ind_auc, "ind_mse" = ind_mse))
}

# A prediction function for GSAM models used in this
# analysis.

gsam_pred <- function(obj, dat) {
  GSAM:::predict.add_mod(obj, new.data = dat,
                         type  = "response")
}
