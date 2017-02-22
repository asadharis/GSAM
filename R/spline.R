
## hello
spline_s <- function(y, y.ord, x.ord, theta, lambda){
  .C("callSS_Fortran", y = as.double(y.ord), 
     x = as.double(x.ord), sy = as.double(theta), 
     lambda = as.double(lambda), n_point = as.integer(length(y)))
}


# A new function to estimae sqrt spline type obj
