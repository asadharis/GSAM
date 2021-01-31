# Main file for running the simulation.

# Load sourcecode from helper files.
source("spamProx.R")
source("wrapperProx.R")
source("genData.R")
source("helpers.R")

run_sim <- function(n = 300, sd.noise = 0.1, seed = 1) {

  #n = 5000; sd.noise = 0.5; seed = 1
  require(microbenchmark)
  require(ggplot2)

  # Generate data.
  dat <- gen_data(seed = seed, n = n, nsim = 1,
                  sd.noise = sd.noise)

  lam1 <- (1e-2)*sqrt(mean(dat$y_mat[,1]^2))
  lam2 <- lam1^2

  # Generate Q matrices for SpAM
  Q_mats <- sapply(c(3,5,10,15,20,30,50),
                   FUN = function(i, x) {
                     ns(x, df = i)
                     }, x = dat$x_mat[,1])

  results <-
  microbenchmark("SSP" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "ssp")
  },
  "spam3" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                         lam1 = lam1, lam2 = lam2, method = "spam",
                         x_matQ = Q_mats[[1]])
  },
  "spam5" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "spam",
                 x_matQ = Q_mats[[2]])
  },"spam10" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "spam",
                 x_matQ = Q_mats[[3]])
  },"spam15" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "spam",
                 x_matQ = Q_mats[[4]])
  },"spam20" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "spam",
                 x_matQ = Q_mats[[5]])
  },"spam30" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "spam",
                 x_matQ = Q_mats[[6]])
  },"spam50" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "spam",
                 x_matQ = Q_mats[[7]])
  },
  "TF0" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "tf",
                 k = 0)
  },
  "TF1" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "tf",
                 k = 1)
  },
  "TF2" = {
    prox_wrapper(r = dat$y_mat[,1], x = dat$x_mat[,1],
                 lam1 = lam1, lam2 = lam2, method = "tf",
                 k = 2)
  }
  )

  return(results)

}


# A function to process the results, generate plots etc.
process_res <- function(sd.noise=0.5, seed = 1) {
  res_n100 <- run_sim(n = 100, sd.noise = sd.noise, seed = seed)
  res_n500 <- run_sim(n = 500, sd.noise = sd.noise, seed = seed)
  res_n1000 <- run_sim(n = 1000, sd.noise = sd.noise, seed = seed)
  res_n5000 <- run_sim(n = 5000, sd.noise = sd.noise, seed = seed)

  require(ggplot2)
  # For the plots we exclude SPAM
  # because it is much faster for obvious reasons.
  p_100 <- myplot(subset(res_n100, expr %in% c("TF0", "TF1", "TF2", "SSP")),
                  title = "Sample size: 100")
  p_500 <- myplot(subset(res_n500, expr %in% c("TF0", "TF1", "TF2", "SSP")),
                  title = "Sample size: 500")
  p_1000 <- myplot(subset(res_n1000, expr %in% c("TF0", "TF1", "TF2", "SSP")),
                  title = "Sample size: 1000")
  p_5000 <- myplot(subset(res_n5000, expr %in% c("TF0", "TF1", "TF2", "SSP")),
                  title = "Sample size: 5000")

  multiplot(p_100, p_500, p_1000, p_5000, cols = 2)

  return(list(res_n100, res_n500, res_n1000, res_n5000))
}

final_results <- process_res()

