plot.res <- function(p=6, scen = 1,vec = c(4,6,7),
                     ord = c(1,2,3)) {
  load(paste0("FinalData/MSEscen",scen,"_p",p,".RData"))
  
  spam.num <- c(3,6,10,20,30,50,80)[vec]
  spam.lab <- paste0("SpAM, M=",spam.num)
  inds.spam <- which(dat$Method==spam.lab[1]| dat$Method==spam.lab[2] |
    dat$Method==spam.lab[3])
  inds.other <- which(dat$Method=="SSP"| dat$Method=="TF, k=0" |
    dat$Method=="TF, k=1" | dat$Method=="TF, k=2" )
  
  dat <- dat[c(inds.spam, inds.other),]
  dat$Method <- factor(dat$Method)
  
  
  lev.low <- levels(dat$Method)==spam.lab[1]
  lev.mod <- levels(dat$Method)==spam.lab[2]
  lev.high <- levels(dat$Method)==spam.lab[3]
  
  levels(dat$Method)[lev.low] <- "SpAM, Low Ord"
  levels(dat$Method)[lev.mod] <- "SpAM, Mod Ord"
  levels(dat$Method)[lev.high] <- "SpAM, High Ord"
  
  dat$Method <- factor(dat$Method, levels(dat$Method)[c(ord,4,5,6,7)])
  
  require(ggplot2)
  g1 <- ggplot(dat) +
    geom_line(aes(x = n, y = MSE, color = Method, linetype = Method), size = 1.2) + 
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=20), legend.text=element_text(size=12),
          legend.title = element_blank())+
    scale_y_log10() + xlab("Sample Size")+
    scale_linetype_manual(values = c(rep(2,3), 3, rep(1,3))) +
    #guides(linetype = guide_legend(keywidth = 4.1)) + 
    guides(color=FALSE, linetype= FALSE) + 
    ggtitle(paste("Scenario", scen))+
    theme(plot.title = element_text(hjust = 0.5,size = 23,face = "bold"))
  
  g1
}
plot.res(6,1,vec = c(2,6,7), ord = c(2,1,3))
plot.res(6,2,vec = c(1,3,6), ord = c(2,1,3))
plot.res(6,3,vec = c(1,2,4), ord = c(2,3,1))
plot.res(6,4,vec = c(1,3,5), ord = c(2,1,3))
plot.res(6,5,vec = c(1,4,6), ord = c(2,1,3))


plot.res(50,1,vec = c(3,5,7), ord = c(1,2,3))
plot.res(50,2,vec = c(1,3,6), ord = c(2,1,3))
plot.res(50,3,vec = c(1,2,4), ord = c(2,3,1))
plot.res(50,4,vec = c(1,3,5), ord = c(2,1,3))
plot.res(50,5,vec = c(1,4,6), ord = c(2,1,3))





plot.object <- function(obj.name = "mod.ssp",
                        scen.num = 4, p = 50, n = 500,
                        title = "SSP", size = 5) {
  
  require(ggplot2)
  source("Models.R")
  require(splines)
  # obj.name = "mod.ssp"
  # scen.num = 4; p = 50; n = 500; title = "SSP"
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
    guides(color=FALSE) + xlab("") + ylab("") + 
    ggtitle(title) +
    theme(axis.text=element_text(size=14))
  # Sample 50 random data sets to plot functions
  set.seed(1)
  samp <- sample(1:100, size = size)
  
  for(i in samp) {
    load(paste0(dir, "/", files[i]))
    fhat <- eval(parse(text = obj.name))$fhat
    fhat <- scale(fhat, scale = FALSE)
    dat.temp <- data.frame("x" = rep(xs, 4), "y" = as.numeric(fhat),
                           "Function" = factor(rep(paste0("Function ", 1:4), each = 1e+3)))
    g1 <- g1 + geom_line(data = dat.temp, aes(x = x, y=y), alpha = 0.4)
  }
  
  g1 + geom_line(size = 1.4)+
    theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"))
}

plot.object("mod.spam3", 4,50,500, "SpAM, M = 3")
plot.object("mod.spam10", 4,50,500, "SpAM, M = 10")
plot.object("mod.spam30", 4,50,500, "SpAM, M = 30")

plot.object("mod.ssp", 4,50,500, "SSP")

plot.object("mod.tf.k0", 4,50,500, "TF, k = 0")
plot.object("mod.tf.k1", 4,50,500, "TF, k = 1")
plot.object("mod.tf.k2", 4,50,500, "TF, k = 2")



plot.object("mod.spam3", 4,6,500, "SpAM, M = 3")
plot.object("mod.spam10", 4,6,500, "SpAM, M = 10")
plot.object("mod.spam30", 4,6,500, "SpAM, M = 30")

plot.object("mod.ssp", 4,6,500, "SSP")

plot.object("mod.tf.k0", 4,6,500, "TF, k = 0")
plot.object("mod.tf.k1", 4,6,500, "TF, k = 1")
plot.object("mod.tf.k2", 4,6,500, "TF, k = 2")

