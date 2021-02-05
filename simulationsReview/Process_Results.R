# In this file we access all the different .RData files for a given simulation setting.
# and aggregate the objects into objects such as plots/tables/data.frames

process.object <- function(scen = 3, p = 6, n = 100) {
  #p = 100; n = 100; scen  = 3

  dir <- paste0("scen", scen, "/p", p, "/n", n)
  files <- paste0(dir,"/", list.files(dir))
  nsim <- length(files)

  mse <- array(NA, dim = c(4,2, nsim))
  i <- 1
  for(file_name in files) {
    #print(i)
    load(file_name)
    mse[,,i] <- as.matrix(fin.mse[,1:2])
    i <- i+1
  }
  mse.mu <- apply(mse, c(1,2), mean)
  mse.se <- apply(mse, c(1,2), function(x){sd(x)/sqrt(length(x))})

  colnames(mse.se) <- colnames(mse.mu) <- c("Paired", "Unpaired")

  labs <- c("SSP", paste0("TF, k=", 0:2) )

  data.frame("MSE" = as.numeric(mse.mu), "MSEse" = as.numeric(mse.se),
             "n" = rep(n, length(mse.mu)), "p" = as.factor(rep(p,length(mse.mu))),
             "Method" = as.factor(rep(labs,2)),
             "Type" = rep(colnames(mse.mu), each = 4) )

}


process.simulation <- function(ns = c(100,300, 500, 700,
                                      900, 1000, 2000, 5000),
                               p=6, scen = 3) {
  #scen.num = 1; p = 6; spam.num = c(2,4,6)
  #ns <- c(100, 300, 500, 800, 1000)
  dat <- data.frame()
  for(i in ns){
    print(i)
    dat <- rbind(dat, process.object(n=i, p=p, scen = scen))
  }
  return(dat)
}


plot.dat <- function(dat, scen = 3, p = 6) {
  require(ggplot2)

    ggplot(data = dat, aes(x = n, y = MSE, color = Type)) +
    geom_line(aes(linetype = Type), size = .8)+
    geom_point(aes(shape = Type), size = 2.5) +
    #geom_errorbar(aes(ymin = MSE-MSEse, ymax = MSE+MSEse))+
    facet_wrap(~Method, scale = "free_y")+
    scale_y_log10()+
    xlab("Sample size (n)") +
      ylab("Mean square error")+
      ggtitle(paste0("Scenario ", scen,"; p = ",p))+
      theme_bw()+
      theme(text = element_text(size = 15)) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))



}

# Plots for p = 6.

dat_p6_scen3 <- process.simulation(ns = c(seq(100, 900, by = 200), 1000),
                                   scen =3, p = 6)
plot.dat(dat_p6_scen3, scen = 3, p = 6)

dat_p6_scen4 <- process.simulation(ns = c(seq(100, 900, by = 200), 1000),
                                   scen =4, p = 6)
plot.dat(dat_p6_scen4, scen = 4, p=6)


dat_p100_scen3 <- process.simulation(ns = c(seq(100,900, by=200), 1000), #seq(1000,2000, by = 200)
                                     scen = 3, p=100)
plot.dat(dat_p100_scen3, scen = 3, p=100)

dat_p100_scen4 <- process.simulation(ns = c(seq(100,900, by=200),1000),
                          scen = 4, p=100)
plot.dat(dat_p100_scen4, scen = 4, p = 100)

