# In this file we access all the different .RData files for a given simulation setting. 
# and aggregate the objects into objects such as plots/tables/data.frames



process.object <- function(obj.name = "mod.spam10",
                           scen.num = 1, p = 6, n = 200,
                           method.spam = FALSE) {
  # obj.name = "mod.spam10";
  # obj.name = "mod.ssp"
  # scen.num = 1; p = 6; n = 200
  # method.spam = 
  
  dir <- paste0("scen", scen.num, "_p", p, "_n", n)
  files <- list.files(dir)
  nsim <- length(files)
  
  ans <- vector("list", nsim)
  for(i in 1:nsim) {
    load(paste0(dir, "/seed", i, ".RData"))
    ans[[i]] <- eval(parse(text = obj.name))
  }
  
  if(method.spam) {
    lam.range <- do.call("rbind",lapply(ans, function(x){range(x$lam)}) )
    lam.rng <- c(min(lam.range[,1]), max(lam.range[,2]))
    lam.seq <- 10^seq(log10(lam.rng[1]), log10(lam.rng[2]), length = 50)
    
    mse <- lapply(ans, function(x){
      approx(x$lam, x$mse.true, xout = lam.seq, rule = 2)$y
    })
    mse <- apply(do.call("cbind", mse), 1, mean)
    lam <- lam.seq/sqrt(n)
    
  } else {
    mse <- apply(do.call("cbind", lapply(ans, function(x){(x$mse.true)})), 1, mean)
    lam <- ans[[1]]$lam
  }
  
  mse.best <- mean(do.call("c", lapply(ans, function(x){(x$mse.true.best)})))
  mse.val <- mean(do.call("c", lapply(ans, function(x){(x$mse.val)})))
  act.val <- mean(do.call("c", lapply(ans, function(x){ x$act.set[x$ind]})))
  list("mse" = mse, "lam" = lam, "mse.best" = mse.best,
       "mse.val" = mse.val, "act.val" = act.val)
}

plot.object <- function(obj.name = "mod.ssp",
                        scen.num = 1, p = 6, n = 200) {
  
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
  
  dat <- data.frame("x" = rep(xs, 4), "y" = as.numeric(f0),
                    "Function" = factor(rep(paste0("Function ", 1:4), each = 1e+3)))
  g1 <- ggplot(dat, aes(x = x, y = y, color = Function)) + 
    geom_line(size = 1.4) + facet_grid(.~Function)+
    guides(color=FALSE)
  # Sample 50 random data sets to plot functions
  set.seed(1)
  samp <- sample(1:500, size = 100)
  
  for(i in samp) {
    load(paste0(dir, "/", files[i]))
    fhat <- eval(parse(text = obj.name))$fhat
    dat.temp <- data.frame("x" = rep(xs, 4), "y" = as.numeric(fhat),
                           "Function" = factor(rep(paste0("Function ", 1:4), each = 1e+3)))
    g1 <- g1 + geom_line(data = dat.temp, aes(x = x, y=y), alpha = 0.1)
  }
   
  g1 + geom_line(size = 1.4)
}

plot.object(n=500)




process.results <- function(scen.num = 1, p = 6, n = 200) {
  #scen.num = 1; p = 6; n = 200
  
  mod.spam3 <- process.object("mod.spam3", scen.num, p, n, TRUE)
  mod.spam6 <- process.object("mod.spam6", scen.num, p, n, TRUE)
  mod.spam10 <- process.object("mod.spam10", scen.num, p, n, TRUE)
  mod.spam20 <- process.object("mod.spam20", scen.num, p, n, TRUE)
  
  mod.tf0 <- process.object("mod.tf.k0", scen.num, p, n, FALSE)
  mod.tf1 <- process.object("mod.tf.k1", scen.num, p, n, FALSE)
  mod.tf2 <- process.object("mod.tf.k2", scen.num, p, n, FALSE)
  
  mod.ssp <- process.object("mod.ssp", scen.num, p, n, FALSE)
  
  mse.all <- c(mod.spam3$mse, mod.spam6$mse, mod.spam10$mse, mod.spam20$mse,
               mod.tf0$mse, mod.tf1$mse, mod.tf2$mse, mod.ssp$mse)
  lams.all <- c(mod.spam3$lam, mod.spam6$lam, mod.spam10$lam, mod.spam20$lam,
                mod.tf0$lam, mod.tf1$lam, mod.tf2$lam, mod.ssp$lam)
  labs <- c(rep(paste0("SpAM, M=", c(3, 6, 10, 20)), each = 50),
            rep(paste0("TF, k=", 0:2), each = 50),
            rep("SSP", 50))
  dat <- data.frame("MSE" = mse.all, "Lambda" = lams.all, "Method" = labs)
    
  require(ggplot2)
  
  
  ggplot(dat) +
    geom_line(aes(x = Lambda, y = MSE, color = Method, linetype = Method), size = 1.2) + 
    scale_x_log10() +
    scale_y_log10() +
    scale_linetype_manual(values = c(rep(3,4), rep(2,1), rep(1,3))) +
    guides(linetype = guide_legend(keywidth = 2.1))
}

# process.results(scen.num = 1, p =6)
# 
# 
# process.results(scen.num = 1, p =100)
# 
# 
# process.results(scen.num = 4, p =6)
# 
# 
# process.results(scen.num = 5, p =100)
# 

# Varying p and function of lambda. 

