
get.data <- function(nvars = 100) {
  load("DataFile.RData")
  
  vars <- apply(dat$x, 2, var)
  ind <- tail(order(vars), nvars)
  
  list("y" = dat$y, "x" = dat$x[,ind])
  
}
