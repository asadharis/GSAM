
run.sim <- function(seed = 1, nvars = 100) {
  require(uniSolve)
  source("spam.R")
  source("helpers.R")
  source("ssp.R")
  source("trendfiltering.R")
  
  dat <- get.data(100)
  rm(riboflavin)
  
  n <- length(dat$y)

  #seed <- 1
  # Obtain index for training set.
  
  # We only use the seed value to split the data into training/test.
  set.seed(seed)
  train <- sample(1:n, 50)
  x.train <- dat$x[train,]
  x.test <- dat$x[-train,]
  y.train <- dat$y[train]
  y.test <- dat$y[-train]
  
  # Obtain the cross-validation folds, we keep the same seed for 
  # this.
  folds <- cut(seq(1,nrow(x.train)), breaks=5,labels=FALSE)
  
  
  # SPAM RESULTS!
  spam3 <- simulation.spam(x.train, y.train, x.test, y.test, 
                              folds = folds, nbasis = 3)
  spam4 <- simulation.spam(x.train, y.train, x.test, y.test, 
                           folds = folds, nbasis = 4)
  spam5 <- simulation.spam(x.train, y.train, x.test, y.test, 
                           folds = folds, nbasis = 5)
  spam6 <- simulation.spam(x.train, y.train, x.test, y.test, 
                           folds = folds, nbasis = 6)
  spam7 <- simulation.spam(x.train, y.train, x.test, y.test, 
                           folds = folds, nbasis = 7)
  
  # SSP RESULTS!
  ssp <- simulation.ssp(x.train, y.train, x.test, y.test, folds)
  
  # Trend filtering results
  tf0 <- simulation.tf(x.train, y.train, x.test, y.test, folds, k=0)
  tf1 <- simulation.tf(x.train, y.train, x.test, y.test, folds, k=1)
  tf2 <- simulation.tf(x.train, y.train, x.test, y.test, folds, k=2)
  
  filename <- paste0("riboflavin_nvars", nvars)

  if(dir.exists(filename)) {
    save(spam3, spam4,spam5, spam6,spam7,  ssp, tf0, tf1, tf2,
         file = paste0(filename, "/seed", seed, ".RData"))
  } else {
    dir.create(filename)
    save(spam3, spam4,spam5, spam6,spam7,  ssp, tf0, tf1, tf2,
         file = paste0(filename, "/seed", seed, ".RData"))
    }
}

args <-  commandArgs(T)
seed <- as.numeric(args[[1]])

nvars <- as.numeric(args[[2]])

run.sim(seed=seed, nvars = nvars)


q(save = "no")
