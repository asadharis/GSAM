
simulation <- function(seed=1, n = 200,
                       num.vars = 100, SNR = 10,
                       scen.num = 1) {
  library(glmgen)
  library(uniSolve)
  source('Generate_Data.R')
  source('Models.R')
  source('spam.R')
  source('ssp.R')
  source('trendfilter.R')

  # seed=1; n = 200;
  # num.vars = 100; SNR = 10;
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
                      SNR = SNR, scenario = scen)

  mod.spam3 <- SimSPAM(dat, p = 3, nlambda = 50)
  mod.spam6<- SimSPAM(dat, p = 6, nlambda = 50)
  mod.spam10 <- SimSPAM(dat, p = 10, nlambda = 50)
  mod.spam20 <- SimSPAM(dat, p = 20, nlambda = 50)

  mod.ssp <- SimSPLINE(dat, lambda.max = 1, lambda.min.ratio = 1e-3)

  mod.tf.k0 <- SimTF(dat, k = 0, lambda.max = 1,
                     lambda.min.ratio = 1e-3)
  mod.tf.k1 <- SimTF(dat, k = 1, lambda.max = 1,
                     lambda.min.ratio = 1e-3)
  mod.tf.k2 <- SimTF(dat, k = 2, lambda.max = 1,
                     lambda.min.ratio = 1e-3)
  # mod.tf.k3 <- SimTF(dat, k = 3, lambda.max = 1,
  #                    lambda.min.ratio = 1e-3)

  # plot(mod.ssp$mse.true,  log = "y")
  # lines(mod.tf.k0$mse.true, col = "blue")
  # lines(mod.tf.k1$mse.true, col = "red")
  # lines(mod.tf.k2$mse.true, col = "green")

  dirname <- paste0("scen", scen.num,"_p", num.vars, "_n", n)
  filename <- paste0(dirname, "/seed", seed, ".RData")

  if(dir.exists(dirname)) {
    save(list = c(paste0("mod.spam", c(3,6,10,20)),
                  "mod.ssp",
                  paste0("mod.tf.k", 0:2)), file = filename)
  } else {
    dir.create(dirname)
    save(list = c(paste0("mod.spam", c(3,6,10,20)),
                  "mod.ssp",
                  paste0("mod.tf.k", 0:2)), file = filename)
  }
}

args <-  commandArgs(T)
seed <- as.numeric(args[[1]])
n <- as.numeric(args[[2]])
num.vars <- as.numeric(args[[3]])
SNR <- as.numeric(args[[4]])
scen.num <- as.numeric(args[[5]])

library(glmgen)
library(uniSolve)
source('Generate_Data.R')
source('Models.R')
source('spam.R')
source('ssp.R')
source('trendfilter.R')

simulation(seed, n, num.vars, SNR, scen.num)
