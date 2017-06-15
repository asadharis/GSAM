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

process.results <- function(scen.num = 1, spam.num = c(2,4,6)) {
  # scen.num = 1; p = 6; spam.num = c(2,4,6)
  
  ans6 <- process.object(scen.num, 6)
  ans100 <- process.object(scen.num, 100)
  ans6$mu <- ans6$mu[c(spam.num, 7:10),]
  ans100$mu <- ans100$mu[c(spam.num, 7:10),]
  
  mse.all <- c(as.numeric(t(ans6$mu)), as.numeric(t(ans100$mu)))
  n.all <- c(rep(as.numeric(colnames(ans6$mu)), nrow(ans6$mu)),
             rep(as.numeric(colnames(ans100$mu)), nrow(ans100$mu)))
  
  spam.vals <- c(3, 6, 10, 20, 30, 50)[spam.num]
  
  labs <- c(rep(paste0("SpAM, M=", spam.vals), each = 20),
            rep("SSP", 20),
            rep(paste0("TF, k=", 0:2), each = 20))
  labs <- c(labs, labs)
  NumbP <- factor(rep(c("p = 6", "p = 100"), each = length(n.all)/2))
  dat <- data.frame("MSE" = mse.all, "n" = n.all, "Method" = labs,
                    "p" = NumbP)
  dat
    
}

plot.results <- function(vec1 = c(1,3,6), vec2 = c(1,3,6),
                         vec3 = c(1,2,3), vec4 = c(1,3,6),
                         vec5 = c(1,3,6)) {
  
  # vec1 = c(1,3,6); vec2 = c(1,3,6);
  # vec3 = c(1,2,3); vec4 = c(1,3,6);
  # vec5 = c(2,3,6)
  
  g1 <- process.results(scen.num = 1, spam.num = vec1)
  
  g2 <- process.results(scen.num = 2, spam.num = vec2)
  
  g3 <- process.results(scen.num = 3, spam.num = vec3)
  
  g4 <- process.results(scen.num = 4, spam.num = vec4)
  
  g5 <- process.results(scen.num = 5, spam.num = vec5)
  
  scen <- rep(paste("Scenario", 1:5), each = nrow(g1))
  
  dat <- cbind(rbind(g1,g2,g3,g4,g5),"scen" = scen)
  
  
  
  dat$Method <- factor(dat$Method, levels(dat$Method)[c(2,8,1,3,4,5,6,7)])
  
  require(ggplot2)
  g1 <- ggplot(subset(dat, p == "p = 6")) +
    geom_line(aes(x = n, y = MSE, color = Method, linetype = Method), size = 1.2) + 
    theme(axis.text=element_text(size=12), 
          axis.title=element_text(size=14), legend.text=element_text(size=12),
          legend.title = element_blank())+
    scale_y_log10() +
    scale_linetype_manual(values = c(rep(3,4), 2, rep(1,3))) +
    #geom_vline(aes(xintercept=vl), data=vline.dat)+
    guides(linetype = guide_legend(keywidth = 4.1)) + 
    facet_wrap(~scen, ncol = 2,scales = "free_y") + ggtitle("Low-dimensional Setting")+
    theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold"))+
    theme(legend.position=c(0.7, 0.14))
  
  g2 <- ggplot(subset(dat, p == "p = 100")) +
    geom_line(aes(x = n, y = MSE, color = Method, linetype = Method), size = 1.2) + 
    theme(axis.text=element_text(size=12), 
          axis.title=element_text(size=14), legend.text=element_text(size=12),
          legend.title = element_blank())+
    scale_y_log10() +
    scale_linetype_manual(values = c(rep(3,4), 2, rep(1,3))) +
    #geom_vline(aes(xintercept=vl), data=vline.dat)+
    guides(linetype = guide_legend(keywidth = 4.1)) + 
    facet_wrap(~scen, ncol = 2,scales = "free_y") + ggtitle("High-dimensional Setting")+
    theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold"))+
    theme(legend.position=c(0.7, 0.14))
  list(g1,g2)
}

check.lam.spp <- function(scen.num=1, p=6) {
  ans <- process.object(scen.num, p)
  myans <- ans$mu[nrow(ans$mu),]
  ms <- seq(50,1000, by = 50)
  plot(ms,myans, ylab = "Lambda", xlab = "n")
  
}

#process.results(scen.num = 1, p =6, spam.num = c(2,4,6))

