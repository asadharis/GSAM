
get.data <- function(nvars = 100) {
  load("DataFile.RData")
  
  vars <- apply(dat[,-1], 2, var)
  ind <- tail(order(vars), nvars)
  
  list("y" = dat[,1], "x" = dat[,ind + 1])
  
}
