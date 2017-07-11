
get.data <- function(nvars = 100) {
  load("data/BreastCancerProcessed.RData")
  # vars <- apply(dat[,-1], 2, var)
  # ind <- tail(order(vars), nvars)

  list("y" = dat[,1], "x" = dat[,-1])

}
