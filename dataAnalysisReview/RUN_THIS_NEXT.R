process.dat <- function(nvars = 100) {
  filename <- paste0("riboflavin_nvars", nvars)
  files <- paste0(filename, "/", list.files(filename))
  nfiles <- length(files)
  err <- matrix(NA, ncol = 7, nrow = nfiles)
  colnames(err) <- c(paste0("spam",c(1,3,6)), "ssp", paste0("tf",0:2))
  sparse <- err
  for(i in 1:nfiles) {
    load(files[i])
    err[i, ] <- c(spam1$err["min"],spam3$err["min"], 
                  spam6$err["min"], ssp$err["min"],
                  tf0$err["min"], tf1$err["min"], tf2$err["min"])
    
    sparse[i, ] <- 1 - c(length(spam1$sparse$min),
                     length(spam3$sparse$min),length(spam6$sparse$min),
                     length(ssp$sparse$min),length(tf0$sparse$min),
                     length(tf1$sparse$min),length(tf2$sparse$min))/nvars
  }
  
  error <- as.numeric(err)
  sparsity <- as.numeric(sparse)
  labs <- c(rep("Lasso", 100), rep("SpAM, M=3",100),
            rep("SpAM, M=6",100), rep("SSP",100),
            rep("TF, k=0",100), rep("TF, k=1",100), 
            rep("TF, k=2",100))
  
  mydat <- data.frame("Error" = error,
             "Sparsity" = sparsity, 
             "Method" = factor(labs))
  require(ggplot2)
  
  plot1 <- ggplot(mydat, aes(x=Method, y=Error, fill = Method)) + 
    theme(text = element_text(size=15)) + 
    ylab("Test Error") + xlab("") + guides(fill = FALSE) +
    theme(axis.text=element_text(size=14, face = "bold")) +
    geom_boxplot()+scale_y_log10()
  
  plot2 <- ggplot(mydat, aes(x=Method, y=Sparsity, fill = Method)) + 
    theme(text = element_text(size=15)) + 
    ylab("Sparsity") + xlab("") + guides(fill = FALSE) +
    theme(axis.text=element_text(size=14, face = "bold")) +
    geom_boxplot()
  
  list("p1" = plot1, "p2" = plot2)
}

make.example.plots <- function(nvars = 100, fi=1) {
  filename <- paste0("riboflavin_nvars", nvars)
  files <- paste0(filename, "/", list.files(filename))
  nfiles <- length(files)
  
  # Find the number of common functions
  # We want to look for cases with 4 common plots
  
  num.common <- numeric(100)
  for(i in 1:100){
    load(files[i])
    a <- intersect(spam1$sparse$min,spam3$sparse$min)
    b <- intersect(ssp$sparse$min,tf2$sparse$min)
    
    num.common[i] <- length(intersect(a,b))
  }
  ind <- which(num.common == 4)[1]
  load(files[ind])
  a <- intersect(spam1$sparse$min,spam3$sparse$min)
  b <- intersect(ssp$sparse$min,tf2$sparse$min)
  fin <- intersect(a,b)
  
  spam1.ind <- which(spam1$sparse$min %in% fin)
  spam3.ind <- which(spam3$sparse$min %in% fin)
  ssp.ind <- which(ssp$sparse$min %in% fin)
  tf2.ind <- which(tf2$sparse$min %in% fin)
  
  spam1.fhat <- build.func(spam1$min.plot$xmat[, spam1.ind],
             spam1$min.plot$yhat[, spam1.ind])
  spam3.fhat <- build.func(spam3$min.plot$xmat[, spam3.ind],
                           spam3$min.plot$yhat[, spam3.ind])
  ssp.fhat <- build.func(ssp$min.plot$xmat[, ssp.ind],
                           ssp$min.plot$yhat[, ssp.ind])
  tf.fhat <- build.func(tf2$min.plot$xmat[, tf2.ind],
                        tf2$min.plot$yhat[, tf2.ind])
  cols <- gg_color_hue(7)
  
  # 1 2 4 7
  par(mar = c(5, 6, 4, 2) + 0.1)
  
  # Build plot 1
  rng <- range(c(spam1.fhat$fhat[,fi],spam3.fhat$fhat[,fi],
           ssp.fhat$fhat[,fi], tf.fhat$fhat[,fi]))
  if(fi==1){
    plot(spam1.fhat$x[,fi], spam1.fhat$fhat[,fi],col = cols[1], 
         type = "l", lwd = 2, xlab = "x", ylab = expression('f'[1]*'(x)'),
         cex.lab = 2, ylim = rng)
  }
  if(fi==2){
    plot(spam1.fhat$x[,fi], spam1.fhat$fhat[,fi],col = cols[1], 
         type = "l", lwd = 2, xlab = "x", ylab = expression('f'[15]*'(x)'),
         cex.lab = 2, ylim = rng)
  }
  if(fi==3){
    plot(spam1.fhat$x[,fi], spam1.fhat$fhat[,fi],col = cols[1], 
         type = "l", lwd = 2, xlab = "x", ylab = expression('f'[50]*'(x)'),
         cex.lab = 2, ylim = rng)
  }
  if(fi==4){
    plot(spam1.fhat$x[,fi], spam1.fhat$fhat[,fi],col = cols[1], 
         type = "l", lwd = 2, xlab = "x", ylab = expression('f'[74]*'(x)'),
         cex.lab = 2, ylim = rng)
  }
  
  
  
  lines(spam3.fhat$x[,fi], spam3.fhat$fhat[,fi], col=cols[2],
        lwd = 2)
  lines(ssp.fhat$x[,fi], ssp.fhat$fhat[,fi], col=cols[4],
        lwd = 2)
  lines(tf.fhat$x[,fi], tf.fhat$fhat[,fi], col=cols[7],
        lwd = 2)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

build.func <- function(xmat, yhat) {
  x1 <- seq(range(xmat[,1])[1],range(xmat[,1])[2] , length = 1e+3)
  x2 <- seq(range(xmat[,2])[1],range(xmat[,2])[2], length = 1e+3)
  x3 <- seq(range(xmat[,3])[1],range(xmat[,3])[2], length = 1e+3)
  x4 <- seq(range(xmat[,4])[1],range(xmat[,4])[2], length = 1e+3)
  
  fhat <- matrix(0, ncol = 4, nrow = 1e+3)
  fhat[,1] <- approx(xmat[,1], yhat[,1], xout = x1, rule = 2)$y
  fhat[,2] <- approx(xmat[,2], yhat[,2], xout = x2, rule = 2)$y
  fhat[,3] <- approx(xmat[,3], yhat[,3], xout = x3, rule = 2)$y
  fhat[,4] <- approx(xmat[,4], yhat[,4], xout = x4, rule = 2)$y
  list("fhat" = fhat, "x" = cbind(x1,x2,x3,x4))
}

plot(0:5,0:5, type = "n")
cols <- gg_color_hue(7)
abline(h = 1, col = cols[1], lwd = 2)
abline(h = 2, col = cols[2], lwd = 2)
abline(h = 3, col = cols[4], lwd = 2)
abline(h = 4, col = cols[7], lwd = 2)
