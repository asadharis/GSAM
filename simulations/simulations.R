
simulation <- function(seed=1, n = 100,
                       num.vars = 6, noise.var = 1,
                       scen.num = 1) {
  library(glmgen)
  library(uniSolve)
  source('Generate_Data.R')
  source('Models.R')
  source('spam.R')
  source('ssp.R')
  source('trendfilter.R')

  # n = 100; seed =1
  # num.vars = 6; noise.var = 1;
  # scen.num <- 1

  if(scen.num == 1){
    scen = scen1
  } else if(scen.num == 2){
    scen = scen2
  } else if(scen.num == 3){
    scen = scen3
  } else if(scen.num == 4){
    scen = scen4
  } else if(scen.num == 5){
    scen = scen5
  }

  dat <- GenerateData(seed = seed, n = n, p = num.vars,
                      noise.var = noise.var, scenario = scen)

  mod.spam3 <- SimSPAM(dat, p = 3, nlambda = 50, lambda.min.ratio = 5e-4)
  mod.spam6<- SimSPAM(dat, p = 6, nlambda = 50, lambda.min.ratio = 5e-4)
  mod.spam10 <- SimSPAM(dat, p = 10, nlambda = 50, lambda.min.ratio = 5e-4)
  mod.spam20 <- SimSPAM(dat, p = 20, nlambda = 50, lambda.min.ratio = 5e-4)
  mod.spam30 <- SimSPAM(dat, p = 30, nlambda = 50, lambda.min.ratio = 5e-4)
  mod.spam50 <- SimSPAM(dat, p = 50, nlambda = 50, lambda.min.ratio = 5e-4)
  mod.spam100 <- SimSPAM(dat, p = 100, nlambda = 50, lambda.min.ratio = 5e-4)

  mod.ssp <- SimSPLINE(dat, lambda.max = 1, lambda.min.ratio = 1e-2,
                       tol = 1e-5, max.iter = 3000)

  mod.tf.k0 <- SimTF(dat, k = 0, lambda.max = 3,
                     lambda.min.ratio = 1e-2, tol = 1e-5, max.iter = 3000)
  mod.tf.k1 <- SimTF(dat, k = 1, lambda.max = 1,
                     lambda.min.ratio = 1e-2, tol = 1e-5, max.iter = 3000)
  mod.tf.k2 <- SimTF(dat, k = 2, lambda.max = 1,
                     lambda.min.ratio = 1e-3, tol = 1e-5, max.iter = 3000,
                     control = trendfilter.control.list(obj_tol = 1e-12, max_iter = 600))

  dirname <- paste0("scen", scen.num,"_p", num.vars,"_n",n)
  filename <- paste0(dirname, "/",seed, ".RData")

  if(dir.exists(dirname)) {
    save(list = c(paste0("mod.spam", c(3, 6, 10, 20, 30, 50, 100)),
                  "mod.ssp",
                  paste0("mod.tf.k", 0:2)), file = filename)
  } else {
    dir.create(dirname)
    save(list = c(paste0("mod.spam", c(3, 6, 10, 20, 30, 50, 100)),
                  "mod.ssp",
                  paste0("mod.tf.k", 0:2)), file = filename)
  }
}

args <-  commandArgs(T)
seed <- as.numeric(args[[1]])
n <- as.numeric(args[[2]])
num.vars <- as.numeric(args[[3]])
noise.var <- as.numeric(args[[4]])
scen.num <- as.numeric(args[[5]])

library(glmgen)
library(uniSolve)
source('Generate_Data.R')
source('Models.R')
source('spam.R')
source('ssp.R')
source('trendfilter.R')

simulation(seed, n, num.vars, noise.var, scen.num)
q(save = "no")
