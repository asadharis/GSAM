# In this file we access all the different .RData files for a given simulation setting.
# and aggregate the objects into objects such as plots/tables/data.frames

process.object <- function(scen.num = 1, p = 6, n = 100,
                           spam.num = c(2,4,6)) {
  scen.num = 1; p = 6; n = 1000;spam.num = c(2,4,6)
  
  dir <- paste0("scen", scen.num, "_p", p, "_n", n)
  files <- paste0(dir,"/", list.files(dir))
  nsim <- length(files)
  spam.vals <- c(3,6,10,20,30,50)[spam.num]
  
  mse <- matrix(0, ncol = nsim, nrow = 7)
  for(i in 1:nsim) {
    load(files[i])
    temp.spam <- mget(paste0("mod.spam",spam.vals))
    mse[,i] <- c(mod.ssp$mse.true.best, mod.tf.k0$mse.true.best,
             mod.tf.k1$mse.true.best, mod.tf.k2$mse.true.best,
             temp.spam[[1]]$mse.true.best, temp.spam[[2]]$mse.true.best,
             temp.spam[[3]]$mse.true.best)
  }
  mse.mu <- apply(mse, 1, mean)
  mse.se <- apply(mse, 1, function(x){sd(x)/sqrt(length(x))})
  
  labs <- c("SSP", paste0("TF, k=", 0:2) , paste0("SpAM, M=", spam.vals))

  data.frame("MSE" = mse.mu, "MSEse" = mse.se,
             "n" = rep(n, 7), "p" = as.factor(rep(p, 7)), 
             "Method" = as.factor(labs))

}

process.simulation <- function(scen.num = 1, p = 6, 
                               spam.num = c(2,4,6)) {
  scen.num = 1; p = 6; spam.num = c(2,4,6)
  dat <- data.frame()
  for(i in seq(100, 1000, by = 100)){
    print(i)
    dat <- rbind(dat, process.object(scen.num, p, n=i, spam.num = spam.num))
  }
  if(dir.exists("FinalData")){
    save(dat, file = paste0("FinalData/scen",scen.num, "_p",p,".RData"))
  }else{
    dir.create("FinalData")
    save(dat, file = paste0("FinalData/scen",scen.num, "_p",p,".RData"))
  }
}

############################################################
############################################################
############################################################
############################################################


plot.object <- function(obj.name = "mod.ssp",
                        scen.num = 1, p = 6, n = 200,
                        title = "SSP") {

  require(ggplot2)
  source("Models.R")
  # obj.name = "mod.ssp"
  # scen.num = 1; p = 6; n = 200
  dir <- paste0("scen", scen.num, "_p", p, "_n", n)
  files <- list.files(dir)
  nsim <- length(files)
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

  xs <- seq(-2.5, 2.5, length = 1e+3)
  f0 <- scen(cbind(xs,xs,xs,xs))
  f0 <- scale(f0, scale = FALSE)

  dat <- data.frame("x" = rep(xs, 4), "y" = as.numeric(f0),
                    "Function" = factor(rep(paste0("Function ", 1:4), each = 1e+3)))
  g1 <- ggplot(dat, aes(x = x, y = y, color = Function)) +
    geom_line(size = 1.4) + facet_grid(.~Function)+
    guides(color=FALSE) + xlab("x") + ylab("f(x)") + ggtitle(title)
  # Sample 50 random data sets to plot functions
  set.seed(1)
  samp <- sample(1:500, size = 10)

  for(i in samp) {
    load(paste0(dir, "/", files[i]))
    fhat <- eval(parse(text = obj.name))$fhat
    fhat <- scale(fhat, scale = FALSE)
    dat.temp <- data.frame("x" = rep(xs, 4), "y" = as.numeric(fhat),
                           "Function" = factor(rep(paste0("Function ", 1:4), each = 1e+3)))
    g1 <- g1 + geom_line(data = dat.temp, aes(x = x, y=y), alpha = 0.2)
  }

  g1 + geom_line(size = 1.4)+
    theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"))
}
