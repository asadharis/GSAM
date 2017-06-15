
get.data <- function(nvars = 100) {
  require(qut)
  data("riboflavin")
  
  vars <- apply(riboflavin$x, 2, var)
  ind <- tail(order(vars), nvars)
  
  list("y" = riboflavin$y, "x" = riboflavin$x[,ind])
  
}
