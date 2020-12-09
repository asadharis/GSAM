# A test file to start.

library(readr)
library(tibble)
library(dplyr)
library(pROC)

dat_full <- read_csv("Breast_GSE70947.csv")
n <- nrow(dat_full)

ind_norm <- which(dat_full$type =="normal")
ind_cancer <- (1:n)[-ind_norm]
cancer_yes <- numeric(n)
cancer_yes[ind_cancer] <- 1

# Let's start with a PCA of the variables.
#
# pca_full <- prcomp(dat_full[,-(1:2)], rank. = 100)
# summary(pca_full)
#
# new_cov <- pca_full$x

# Generate a small data with the first 10 principal components.
### NOT TESTING WITH THIS, PCA components too predictive.
# mydat <- tibble("y" = as.factor(cancer_yes)) %>%
#   bind_cols(as_tibble(new_cov[,1:10]))

# Generate a small data picking the first 10 genes.
mydat <- tibble("y" = as.factor(cancer_yes)) %>%
  bind_cols(dat_full[,3:102])

# plot(new_cov[ind_norm,1], new_cov[ind_norm,2], col = "blue",
#      xlim = range(new_cov[,1]), ylim = range(new_cov[,2]))
# points(new_cov[ind_cancer,1], new_cov[ind_cancer,2], col = "red")

# Let's create a simple training, testing and validation set.
set.seed(1)
sam_set <- sample(1:3, size = n, replace=TRUE, prob = c(0.60,0.2,0.2))
dat_train <- mydat[sam_set==1,]
dat_test <- mydat[sam_set==2,]
dat_val <- mydat[sam_set==3,]

# Fit a simple GLM
mod_glm <- glm(y~., data = dat_train,
               family = binomial())
test_phat <- predict(mod_glm, type = "response", newdata = dat_test)
val_phat <- predict(mod_glm, type = "response", newdata = dat_val)

roc_train <- roc(response = dat_train$y, predictor = mod_glm$fitted.values)
roc_test <- roc(response = dat_test$y, predictor = test_phat)
roc_val <- roc(response = dat_val$y, predictor = val_phat)

cat("Train AUC: ", round(roc_train$auc*100,1), "%\n",sep = "")
cat("Test AUC: ", round(roc_test$auc*100,1), "%\n",sep = "")
cat("Val AUC: ", round(roc_val$auc*100,1), "%\n",sep = "")

###########################################################
###########################################################
###########################################################

# A helper function for running analysis from models.
run.test <- function(obj, dat_train, dat_test, dat_val,
                     pred.fun = predict.spamLL) {
  train_phat <- pred.fun(obj, as.matrix(dat_train[,-1]))
  test_phat <- pred.fun(obj, as.matrix(dat_test[,-1]))
  val_phat <- pred.fun(obj, as.matrix(dat_val[,-1]))

  auc_train <- auc_test <- auc_val <- numeric(ncol(train_phat))
  for(i in 1:ncol(train_phat) ){
    #print(i)
    roc_train <- roc(response = dat_train$y, predictor = train_phat[,i], quiet = TRUE)
    roc_test <- roc(response = dat_test$y, predictor = test_phat[,i], quiet = TRUE)
    roc_val <- roc(response = dat_val$y, predictor = val_phat[,i], quiet = TRUE)
    auc_train[i] <- roc_train$auc
    auc_test[i] <- roc_test$auc
    auc_val[i] <- roc_val$auc
  }

  ind <- which.max(auc_test)
  cat("Train AUC: ", round(auc_train[ind]*100,1), "%\n",sep = "")
  cat("Test AUC: ", round(auc_test[ind]*100,1), "%\n",sep = "")
  cat("Val AUC: ", round(auc_val[ind]*100,1), "%\n",sep = "")
  return(c("train" = auc_train[ind], "test" = auc_test[ind],
           "val" = auc_val[ind]))
}

###########################################################
###########################################################
###########################################################
# Lasso now.

library(glmnet)
mod_lasso <- glmnet(x = as.matrix(dat_train[,-1]),y =  dat_train$y,
                    family = "binomial",
                    nlambda = 50, lambda.min.ratio = 1e-5)
pred.lasso <- function(obj, dat) {
  predict(obj, newx = dat, type = "response")
}
lasso_res <- run.test(obj = mod_lasso,dat_train =  dat_train,dat_test =  dat_test,
         dat_val = dat_val, pred.fun = pred.lasso)


###########################################################
###########################################################
###########################################################
# Now SPAM, sparse additive modeling
library(SAM)
#?samLL

mod_spam2 <- samLL(X = dat_train[,-1],
                   y = as.numeric(as.character(dat_train$y)), p=2,
                   nlambda = 50, lambda.min.ratio = 1e-2)
mod_spam4 <- samLL(X = dat_train[,-1],
                   y = as.numeric(as.character(dat_train$y)), p=4,
                   nlambda = 50, lambda.min.ratio = 1e-2)
mod_spam6 <- samLL(X = dat_train[,-1],
                   y = as.numeric(as.character(dat_train$y)), p=6,
                   nlambda = 50, lambda.min.ratio = 1e-2)
spam2_res <- run.test(mod_spam2, dat_train, dat_test, dat_val,
                      pred.fun = predict.spamLL)
spam4_res <- run.test(mod_spam4, dat_train, dat_test, dat_val,
                      pred.fun = predict.spamLL)
spam6_res <- run.test(mod_spam6, dat_train, dat_test, dat_val,
                      pred.fun = predict.spamLL)

######################################################################
######################################################################
######################################################################

# Now for GSAM with Sobolev penalty

library(GSAM)
mod_ssp <- fit.additive(y=dat_train$y , x= dat_train[, -1],
                         family="binomial", method = "sobolev",
                         lambda.max = 3,
                         lambda.min.ratio = 1e-3, max.iter = 100)

gsam_pred <- function(obj, dat) {
  GSAM:::predict.add_mod(obj, new.data = dat,
                         type  = "response")
}

ssp_res <- run.test(mod_ssp, dat_train, dat_test, dat_val, pred.fun = gsam_pred)


######################################################################
######################################################################
######################################################################

# Finally, we have the Trend filtering penalties.


mod_tf0 <- fit.additive(y=dat_train$y , x= dat_train[, -1],
                        family="binomial", method = "tf", k = 0,
                        lambda.max = 3,
                        lambda.min.ratio = 1e-3, max.iter = 100)

mod_tf1 <- fit.additive(y=dat_train$y , x= dat_train[, -1],
                        family="binomial", method = "tf", k = 1,
                        lambda.max = 3,
                        lambda.min.ratio = 1e-3, max.iter = 100)
mod_tf2 <- fit.additive(y=dat_train$y , x= dat_train[, -1],
                        family="binomial", method = "tf", k = 2,
                        lambda.max = 3,
                        lambda.min.ratio = 1e-3, max.iter = 100)
tf0_res <- run.test(mod_tf0, dat_train, dat_test, dat_val,
                    pred.fun = gsam_pred)
tf1_res <- run.test(mod_tf1, dat_train, dat_test, dat_val,
                    pred.fun = gsam_pred)
tf2_res <- run.test(mod_tf2, dat_train, dat_test, dat_val,
                    pred.fun = gsam_pred)



########################################################
########################################################
########################################################

# Testing parallelization of method
max.lambda <- 1
lambda.min.ratio <- 1e-3

mod_serial <- fit.additive(y=dat_train$y, x=dat_train[,-1],
                           family="binomial", method = "sobolev",
                           lambda.max = max.lambda,
                           lambda.min.ratio = lambda.min.ratio,
                           max.iter = 100, parallel = FALSE)

mod_para <- fit.additive(y=dat_train$y, x=dat_train[,-1],
                        family="binomial", method = "sobolev",
                        lambda.max = max.lambda,
                        lambda.min.ratio = lambda.min.ratio,
                        max.iter = 100, parallel = TRUE)
