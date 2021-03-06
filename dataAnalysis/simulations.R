
run.sim <- function(seed = 1) {
  require(uniSolve)
  source("spam.R")
  source("lasso.R")

  source("ssp.R")
  source("trendfiltering.R")

  load("Data/BostonHousing.RData")
  n <- length(dat$y)

  #
  seed = 1
  # We only use the seed value to split the data into training/test.
  set.seed(seed)
  n <- length(dat$y)
  train <- sample(1:n, floor(n*.75))
  x.train <- as.matrix(dat$x[train,])
  x.test <- as.matrix(dat$x[-train,])
  y.train <- dat$y[train]
  y.test <- dat$y[-train]

  # # Obtain the cross-validation folds, we keep the same seed for
  # # this.
  folds <- cut(seq(1,nrow(x.train)), breaks=5,labels=FALSE)

  # Lasso Results First
  lasso <- simulation.lasso(x.train, y.train, x.test, y.test, folds,
                            lambda.min.ratio = 1e-3)

  # SPAM RESULTS!


  spam2 <- simulation.spam(x.train, y.train, x.test, y.test,
                           folds = folds, nbasis = 2)
  spam3 <- simulation.spam(x.train, y.train, x.test, y.test,
                              folds = folds, nbasis = 3)
  # spam4 <- simulation.spam(x.train, y.train, x.test, y.test,
  #                          folds = folds, nbasis = 4)
  # spam5 <- simulation.spam(x.train, y.train, x.test, y.test,
  #                          folds = folds, nbasis = 5)
  # spam8 <- simulation.spam(x.train, y.train, x.test, y.test,
  #                          folds = folds, nbasis = 8)
  # spam10 <- simulation.spam(x.train, y.train, x.test, y.test,
  #                          folds = folds, nbasis = 10)


  # SSP RESULTS!
  ssp <- simulation.ssp(x.train, y.train, x.test, y.test, folds,
                        max.lambda = 1, lam.min.ratio = 1e-2)

  # Trend filtering results
  tf0 <- simulation.tf(x.train, y.train, x.test, y.test, folds, k=0,
                       lambda.min.ratio = 1e-1,lambda.max = 1)

  tf1 <- simulation.tf(x.train, y.train, x.test, y.test, folds, k=1,
                       lambda.min.ratio = 1e-4,lambda.max = 1)

  tf2 <- simulation.tf(x.train, y.train, x.test, y.test, folds, k=2,
                       lambda.min.ratio = 1e-4,lambda.max = 1)


  filename <- paste0("ERdata")

  if(dir.exists(filename)) {
    save(lasso, spam2, spam3, #spam4,spam5, spam6,
         ssp, tf0, tf1, tf2,
         file = paste0(filename, "/seed", seed, ".RData"))
  } else {
    dir.create(filename)
    save(spam1, spam2, spam3, #spam4,spam5, spam6,
         ssp, tf0, tf1, tf2,
         file = paste0(filename, "/seed", seed, ".RData"))
    }
}

args <-  commandArgs(T)
seed <- as.numeric(args[[1]])

run.sim(seed=seed)


q(save = "no")
