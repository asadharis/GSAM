tf_dp <- function(y, y.ord, lambda, theta) {
  .C("tf_dp_R",
        n = as.integer(length(y)),
        y = as.double(y.ord),
        lam1 = as.double(lambda),
        beta = as.double(theta))
}


tf_test <- function(y, y.ord, lambda) {
  n <- length(y)
  .C("tf_wrap",
     n = as.integer(n),
     y = as.double(y.ord),
     lam = as.double(lambda),
     beta = as.double(rep(0, n)))
}

