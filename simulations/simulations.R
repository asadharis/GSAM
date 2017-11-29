
simulation <- function(n = 100,
                       num.vars = 6, noise.var = 1,
                       scen.num = 1) {
  library(glmgen)
  library(uniSolve)
  source('Generate_Data.R')
  source('Models.R')
  source('spam.R')
  source('ssp.R')
  source('trendfilter.R')

  # n = 100;
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
  mse.vals <- matrix(NA, nrow = 500, ncol = 11)
  colnames(mse.vals) <- c(paste0("spam", c(3,6,10,20,30,50,100)),
                          "ssp", paste0("tf",0:2) )
  sp3 <- sp6 <- sp10 <- sp20 <- sp30 <- sp50 <- sp100 <-
    ssp <- tf0 <- tf1 <- tf2 <- array(NA, dim = c(1000,4,500))
  for(seed in 1:500){
    print(seed)
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
    mse.vals[seed,] <- c(mod.spam3$mse.true.best,mod.spam6$mse.true.best,
                         mod.spam10$mse.true.best,mod.spam20$mse.true.best,
                         mod.spam30$mse.true.best,mod.spam50$mse.true.best,
                         mod.spam100$mse.true.best, mod.ssp$mse.true.best,
                         mod.tf.k0$mse.true.best,mod.tf.k1$mse.true.best,
                         mod.tf.k2$mse.true.best)
    sp3[,,seed] <- mod.spam3$fhat
    sp6[,,seed] <- mod.spam3$fhat
    sp10[,,seed] <- mod.spam3$fhat
    sp20[,,seed] <- mod.spam3$fhat
    sp30[,,seed] <- mod.spam3$fhat
    sp50[,,seed] <- mod.spam3$fhat
    sp100[,,seed] <- mod.spam3$fhat
    ssp[,,seed] <- mod.ssp$fhat
    tf0[,,seed] <- mod.tf.k0$fhat
    tf1[,,seed] <- mod.tf.k1$fhat
    tf2[,,seed] <- mod.tf.k2$fhat
  }

  dirname <- paste0("scen", scen.num,"_p", num.vars)
  filename <- paste0(dirname, "/n", n, ".RData")

  if(dir.exists(dirname)) {
    save(list = c(paste0("sp", c(3, 6, 10, 20, 30, 50, 100)),
                  "ssp",
                  paste0("tf", 0:2), mse.vals), file = filename)
  } else {
    dir.create(dirname)
    save(list = c(paste0("sp", c(3, 6, 10, 20, 30, 50, 100)),
                  "ssp",
                  paste0("tf", 0:2), mse.vals), file = filename)
  }
}

args <-  commandArgs(T)
n <- as.numeric(args[[1]])
num.vars <- as.numeric(args[[2]])
noise.var <- as.numeric(args[[3]])
scen.num <- as.numeric(args[[4]])

library(glmgen)
library(uniSolve)
source('Generate_Data.R')
source('Models.R')
source('spam.R')
source('ssp.R')
source('trendfilter.R')

simulation(seed, n, num.vars, noise.var, scen.num)
q(save = "no")
