# In this file we access all the different .RData files for a given simulation setting. 
# and aggregate the objects into objects such as plots/tables/data.frames

process.object <- function(scen.num = 1, p = 6) {
  # scen.num = 1; p = 6;
  dir <- paste0("sim2scen", scen.num, "_p", p)
  load(paste0(dir, "/seed", 1, ".RData"))
  files <- list.files(dir)
  nsim <- length(files)
  ans <- array(NaN, dim = c(dim(ans.mat), nsim))
  for(i in 1:nsim) {
    load(paste0(dir, "/", files[i]))
    ans[,,i] <- ans.mat
  }
  dimnames(ans) <- list(rownames(ans.mat),
                        colnames(ans.mat),
                        rep("", nsim))
  
  mu.ans <- apply(ans, c(1,2), mean)
  mu.se <- apply(ans, c(1,2), function(x){
    sqrt(var(x)/length(x))
  })
  list("mu" = mu.ans, "se" = mu.se)
}

process.results <- function(scen.num = 1, p = 6, spam.num = c(2,4,6)) {
  # scen.num = 1; p = 6; spam.num = c(2,4,6)
  
  ans <- process.object(scen.num, p)
  ans$mu <- ans$mu[c(spam.num, 7:10),]
  ans$se <- ans$se[c(spam.num, 7:10),]
  
  mse.all <- as.numeric(t(ans$mu))
  n.all <- rep(as.numeric(colnames(ans$mu)), nrow(ans$mu))
  
  spam.vals <- c(3, 6, 10, 20, 30, 50)[spam.num]
  
  labs <- c(rep(paste0("SpAM, M=", spam.vals), each = 20),
            rep("SSP", 20),
            rep(paste0("TF, k=", 0:2), each = 20))
  dat <- data.frame("MSE" = mse.all, "n" = n.all, "Method" = labs)
  
  require(ggplot2)
  ggplot(dat) +
    geom_line(aes(x = n, y = MSE, color = Method, linetype = Method), size = 1.2) + 
    #scale_x_log10() +
    scale_y_log10() +
    scale_linetype_manual(values = c(rep(3,length(spam.num)), 2, rep(1,3))) +
    guides(linetype = guide_legend(keywidth = 2.1))
}

check.lam.spp <- function(scen.num=1, p=6) {
  ans <- process.object(scen.num, p)
  myans <- ans$mu[nrow(ans$mu),]
  ms <- seq(50,1000, by = 50)
  plot(ms,myans, ylab = "Lambda", xlab = "n")
  
}

process.results(scen.num = 1, p =6, spam.num = c(2,4,6))

