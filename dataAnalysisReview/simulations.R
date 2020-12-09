# This file contains the basic simulation function.
# The overall template for a simulation will be as follows:
#   1. With a data split into train, test, and validation.
#   2. For each method:
#       a. Run method on training set for seq of lambda values.
#       b. Fit to test set for all lambda values.
#       c. Calculate test_AUC and test_MSE.
#       d. Find val_AUC and val_MSE for max/min test value.
#       e. Report a warning if test_AUC/MSE at extreme end.
#   3. Aggregate results for each method.
#   4. Save results to DataFile

run.sim <- function(seed = 1, name = "Breast_GSE70947_small",
                    ncores = 8) {
  require(GSAM)
  require(glmnet)
  require(SAM)
  require(glmgen)
  require(parallel)
  require(doParallel)

  source("DataProcessing.R")
  source("helpers.R")
  source("lasso.R")
  source("spam.R")
  source("ssp.R")
  source("trendfiltering.R")

  dat <- process.dat(name = name)

  # Use seed value to split the data into training/test/val.
  set.seed(seed)
  n <- length(dat$y)

  sam_set <- sample(1:3, size = n, replace=TRUE, prob = c(0.60,0.2,0.2))
  dat_train <- dat[sam_set==1,]
  dat_test <- dat[sam_set==2,]
  dat_val <- dat[sam_set==3,]

  # Check if our training set has both labels
  if(length(table(dat_train$y)) != 2){
    stop("Training data does not have both labels.")
  }

  # Lasso Results First -------------------------------
  lasso <- simulation.lasso(dat_train, dat_test, dat_val,
                            lambda.min.ratio = 1e-3)

  cat("Finished lasso\n")

  # SPAM RESULTS -------------------------------
  spam2 <- simulation.spam(dat_train, dat_test, dat_val, nbasis = 2,
                           lambda.min.ratio = 1e-3)
  cat("Finished SPAM2\n")
  spam3 <- simulation.spam(dat_train, dat_test, dat_val, nbasis = 3,
                           lambda.min.ratio = 1e-3)
  cat("Finished SPAM3\n")
  spam4 <- simulation.spam(dat_train, dat_test, dat_val, nbasis = 4,
                           lambda.min.ratio = 1e-3)
  # SPAM larger p -------------------------------------
  cat("Finished SPAM4\n")
  spam5 <- simulation.spam(dat_train, dat_test, dat_val, nbasis = 5,
                           lambda.min.ratio = 1e-3)
  cat("Finished SPAM5\n")
  spam10 <- simulation.spam(dat_train, dat_test, dat_val, nbasis = 10,
                           lambda.min.ratio = 1e-3)

  cat("Finished SPAM10\n")

  # SSP RESULTS!-------------------------------
  ssp <- simulation.ssp(dat_train, dat_test, dat_val,
                        max.lambda = 1, lambda.min.ratio = 1e-3,
                        parallel = TRUE, ncores = ncores)

  cat("Finished SSP\n")
  # Trend filtering results -------------------------------
  tf0 <- simulation.tf(dat_train, dat_test, dat_val, k=0,
                       lambda.min.ratio = 1e-3, max.lambda = 1,
                       parallel = TRUE, ncores = ncores)

  cat("Finished TF0\n")
  # TF 1 -------------------------------
  tf1 <- simulation.tf(dat_train, dat_test, dat_val, k=1,
                       lambda.min.ratio = 1e-3, max.lambda = 1,
                       parallel = TRUE, ncores = ncores)
  cat("Finished TF1\n")
  # TF 2 -------------------------------
  tf2 <- simulation.tf(dat_train, dat_test, dat_val, k=2,
                       lambda.min.ratio = 1e-3, max.lambda = 1,
                       parallel = TRUE, ncores = ncores)
  cat("Finished TF2\n")
  # Collect the results
  method_name <- c("lasso",
                paste0("spam", c(2,3,4,5,10)),
                "ssp",
                paste0("tf", 0:2))

  results <- rbind(lasso, spam2, spam3, spam4, spam5, spam10,
                   ssp, tf0, tf1, tf2)
  rownames(results) <- NULL
  results <- cbind(method_name, results)

  filename <- name

  if(!dir.exists(filename)) {
    dir.create(filename)
  }

  save(results,
       file = paste0(filename, "/seed", seed, ".RData"))

}


# MAIN ############################

 args <-  commandArgs(T)
 seed <- as.numeric(args[[1]])
 name <- as.character(args[[2]])
 ncores <- as.numeric(args[[3]])

#seed <- 1
run.sim(seed=seed, name = name, ncores = ncores)

q(save = "no")
