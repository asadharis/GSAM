tf_dp <- function(y, y.ord, lambda, theta) {
  .C("tf_dp_R",
        n = as.integer(length(y)),
        y = as.double(y.ord),
        lam1 = as.double(lambda),
        beta = as.double(theta))
}
