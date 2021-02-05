# Processing objects so we can draw box plots
process.object2 <- function(scen = 3, p = 6, n = 100) {
  #p = 6; n = 100; scen  = 3

  dir <- paste0("scen", scen, "/p", p, "/n", n)
  files <- paste0(dir,"/", list.files(dir))
  nsim <- length(files)

  mse <- data.frame()
  for(i in 1:nsim) {
    load(files[i])
    mse <- rbind(mse, data.frame("MSE" = c(fin.mse[,1],fin.mse[,2]),
                                 "Method" = rep(fin.mse[,3],2),
                                 "Type" = rep(c("Paired", "Unpaired"),
                                              each = 4)))
  }

  mse$n <- n
  mse$p <- p
  mse
}


process.simulation2 <- function(ns = c(100,300, 500, 700,
                                      900, 1000, 2000, 5000),
                               p=6, scen = 3) {
  #scen.num = 1; p = 6; spam.num = c(2,4,6)
  #ns <- c(100, 300, 500, 800, 1000)
  dat <- data.frame()
  for(i in ns){
    print(i)
    dat <- rbind(dat, process.object2(n=i, p=p, scen = scen))
  }
  return(dat)
}


plot.dat2 <- function(dat) {
  require(ggplot2)

  ggplot(data = dat, aes(x = as.factor(n), y = MSE, color = Type)) +
    geom_boxplot()+
    facet_wrap(~Method,scales = "free_y") +
    scale_y_log10()

}
dat <- process.simulation2(ns = c(100, 300, 500, 700, 900,1000, 2000,5000),
                          scen = 4)
plot.dat2(dat)
