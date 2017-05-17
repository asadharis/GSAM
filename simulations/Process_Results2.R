# In this file we access all the different .RData files for a given simulation setting. 
# and aggregate the objects into objects such as plots/tables/data.frames

process.object <- function(scen.num = 1, p = 6) {
  # obj.name = "mod.spam10";
  # obj.name = "mod.ssp"
  # scen.num = 1; p = 6; n = 200
  
  dir <- paste0("sim2scen", scen.num, "_p", p, "_n200")
  ans <- array(NaN, dim = c(8, 20, 100))
  for(i in 1:100) {
    load(paste0(dir, "/seed", i, ".RData"))
    ans[,,i] <- ans.mat
  }
  dimnames(ans) <- list(rownames(ans.mat), 
                       colnames(ans.mat),
                       rep("", 100))
  
  mu.ans <- apply(ans, c(1,2), mean)
  mu.se <- apply(ans, c(1,2), function(x){
    sqrt(var(x)/length(x))
    })
  list("mu" = mu.ans, "se" = mu.se)
}

process.results <- function(scen.num = 1, p = 6, n = 200) {
  # scen.num = 1; p = 6
  # 
  
  ans <- process.object(scen.num, p)
  
  
  mse.all <- as.numeric(t(ans$mu))
  n.all <- rep(as.numeric(colnames(ans$mu)), 8)
  labs <- c(rep(paste0("SpAM, M=", c(3, 6, 10, 20)), each = 20),
            rep("SSP", 20),
            rep(paste0("TF, k=", 0:2), each = 20))
  dat <- data.frame("MSE" = mse.all, "n" = n.all, "Method" = labs)
  
  require(ggplot2)
  
  
  ggplot(dat) +
    geom_line(aes(x = n, y = MSE, color = Method, linetype = Method), size = 1.2) + 
    #scale_x_log10() +
    scale_y_log10() +
    scale_linetype_manual(values = c(rep(3,4), rep(2,1), rep(1,3))) +
    guides(linetype = guide_legend(keywidth = 2.1))
}

process.results(scen.num = 3, p =6)

