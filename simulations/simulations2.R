
simulation2 <- function(seed=1,
                       num.vars = 100, SNR = 10,
                       scen.num = 1) {
  library(glmgen)
  library(uniSolve)
  source('Generate_Data.R')
  source('Models.R')
  source('spam.R')
  source('ssp.R')
  source('trendfilter.R')

  # seed=1; n.seq = c(50, 100, 200, 500);
  # num.vars = 6; SNR = 10;
  # scen.num <- 5

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

  n.seq <- seq(50, 1000, by = 50)

  ans.mat <- sapply(n.seq, function(n){
    print(n)
    dat <- GenerateData(seed = seed, n = n, p = num.vars,
                        SNR = SNR, scenario = scen)

    mod.spam3 <- SimSPAM(dat, p = 3, nlambda = 50, lambda.min.ratio = 5e-4)$mse.true.best
    mod.spam6<- SimSPAM(dat, p = 6, nlambda = 50, lambda.min.ratio = 5e-4)$mse.true.best
    mod.spam10 <- SimSPAM(dat, p = 10, nlambda = 50, lambda.min.ratio = 5e-4)$mse.true.best
    mod.spam20 <- SimSPAM(dat, p = 20, nlambda = 50, lambda.min.ratio = 5e-4)$mse.true.best

    mod.ssp <- SimSPLINE(dat, lambda.max = 1, lambda.min.ratio = 1e-3)$mse.true.best

    mod.tf.k0 <- SimTF(dat, k = 0, lambda.max = 1,
                       lambda.min.ratio = 1e-3)$mse.true.best
    mod.tf.k1 <- SimTF(dat, k = 1, lambda.max = 1,
                       lambda.min.ratio = 1e-3)$mse.true.best
    mod.tf.k2 <- SimTF(dat, k = 2, lambda.max = 1,
                       lambda.min.ratio = 1e-3)$mse.true.best
    c("spam3" = mod.spam3, "spam6" = mod.spam6,
      "spam10" = mod.spam10, "spam20" = mod.spam20, "ssp" = mod.ssp,
     "tf0" = mod.tf.k0, "tf1" = mod.tf.k1, "tf2" = mod.tf.k2)

  })

  colnames(ans.mat) <- n.seq

  dirname <- paste0("sim2scen", scen.num,"_p", num.vars, "_n", n)
  filename <- paste0(dirname, "/seed", seed, ".RData")

  if(dir.exists(dirname)) {
    save(ans.mat, file = filename)
  } else {
    dir.create(dirname)
    save(ans.mat, file = filename)
  }
}

args <-  commandArgs(T)
seed <- as.numeric(args[[1]])
num.vars <- as.numeric(args[[2]])
SNR <- as.numeric(args[[3]])
scen.num <- as.numeric(args[[4]])

library(glmgen)
library(uniSolve)
source('Generate_Data.R')
source('Models.R')
source('spam.R')
source('ssp.R')
source('trendfilter.R')

simulation2(seed, num.vars, SNR, scen.num)
