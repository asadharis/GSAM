
simulation <- function(seed=1, n = 100,
                       num.vars = 100, noise.var = 1,
                       scen.num = 3) {
  library(glmgen)
  library(GSAM)
  source('Generate_Data.R')
  source('Models.R')
  source('spam.R')
  source('ssp.R')
  source('trendfilter.R')

  # n = 500; seed =1
  # num.vars = 6; noise.var = 1;
  # scen.num <- 3

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

  mod.ssp <- SimSPLINE(dat, lambda.max = 1, lambda.min.ratio = 1e-2,
                       tol = 1e-4, max.iter = 300)

  mod.tf.k0 <- SimTF(dat, k = 0, lambda.max = 2,
                     lambda.min.ratio = 1e-2, tol = 1e-4, max.iter = 300)
  mod.tf.k1 <- SimTF(dat, k = 1, lambda.max = 1,
                     lambda.min.ratio = 1e-3, tol = 1e-4, max.iter = 300)
  mod.tf.k2 <- SimTF(dat, k = 2, lambda.max = 0.1,
                       lambda.min.ratio = 1e-3, tol = 1e-4, max.iter = 300)

  fin.mse <- data.frame(rbind(mod.ssp, mod.tf.k0, mod.tf.k1, mod.tf.k2))
  fin.mse$method <- c("SSP", "TF0", "TF1", "TF2")
  row.names(fin.mse) <- NULL

  dirname <- paste0("p", num.vars,"_n",n)
  filename <- paste0(dirname, "/",seed, ".RData")

  if(dir.exists(dirname)) {
    save(list = c("fin.mse"), file = filename)
  } else {
    dir.create(dirname)
    save(list = c("fin.mse"), file = filename)
  }
}

args <-  commandArgs(T)
seed <- as.numeric(args[[1]])
n <- as.numeric(args[[2]])
num.vars <- as.numeric(args[[3]])
noise.var <- as.numeric(args[[4]])
scen.num <- as.numeric(args[[5]])

library(glmgen)
library(GSAM)
source('Generate_Data.R')
source('Models.R')
source('spam.R')
source('ssp.R')
source('trendfilter.R')

simulation(seed, n, num.vars, noise.var, scen.num)

q(save = "no")
