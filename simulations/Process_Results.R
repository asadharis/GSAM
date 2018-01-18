# In this file we access all the different .RData files for a given simulation setting.
# and aggregate the objects into objects such as plots/tables/data.frames

process.object <- function(scen.num = 1, 
                           p = 6, n = 100) {
  #scen.num = 1; p = 6; n = 800;spam.num = c(2,4,7)
  
  dir <- paste0("scen", scen.num, "_p", p, "_n", n)
  files <- paste0(dir,"/", list.files(dir))
  nsim <- length(files)
  spam.vals <- c(3, 6, 10, 20, 30, 50, 80)#[spam.num]
  
  mse <- ind <- matrix(0, ncol = nsim, nrow = 11)
  for(i in 1:nsim) {
    load(files[i])
    mse[,i] <- c(mod.ssp$mse.true.best, mod.tf.k0$mse.true.best,
             mod.tf.k1$mse.true.best, mod.tf.k2$mse.true.best,
             mod.spam3$mse.true.best, mod.spam6$mse.true.best,
             mod.spam10$mse.true.best,
             mod.spam20$mse.true.best,
             mod.spam30$mse.true.best,mod.spam50$mse.true.best,
             mod.spam80$mse.true.best)
    ind[,i] <- c(mod.ssp$ind, mod.tf.k0$ind,
                 mod.tf.k1$ind, mod.tf.k2$ind,
                 mod.spam3$ind, mod.spam6$ind,
                 mod.spam10$ind,
                 mod.spam20$ind,
                 mod.spam30$ind,mod.spam50$ind,
                 mod.spam80$ind)
  }
  mse.mu <- apply(mse, 1, mean)
  ind.mu <- apply(ind, 1, mean)
  mse.se <- apply(mse, 1, function(x){sd(x)/sqrt(length(x))})
  
  labs <- c("SSP", paste0("TF, k=", 0:2) , paste0("SpAM, M=", spam.vals))

  data.frame("MSE" = mse.mu, "MSEse" = mse.se,
             "n" = rep(n, nrow(mse)), "p" = as.factor(rep(p,nrow(mse))), 
             "Method" = as.factor(labs),"ind"= ind.mu)

}

process.simulation <- function(scen.num = 1, p = 6) {
  #scen.num = 1; p = 6; spam.num = c(2,4,6)
  dat <- data.frame()
  for(i in seq(100, 800, by = 100)){
    print(i)
    dat <- rbind(dat, process.object(scen.num, p, n=i))
  }
  if(dir.exists("FinalData")){
    save(dat, file = paste0("FinalData/MSEscen",scen.num, "_p",p,".RData"))
  }else{
    dir.create("FinalData")
    save(dat, file = paste0("FinalData/MSEscen",scen.num, "_p",p,".RData"))
  }
}

############################################################

process.simulation(1,6)
process.simulation(2,6)
process.simulation(3,6)
process.simulation(4,6)
process.simulation(5,6)
# 
process.simulation(1,50)
process.simulation(2,50)
process.simulation(3,50)
process.simulation(4,50)
process.simulation(5,50)
