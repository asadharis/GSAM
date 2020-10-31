
# This function fits the TF penalty and returns a single
# value: the value of the oracle MSE.
oracle.TF <- function(dat,lambda.max = 1, lambda.min.ratio = 1e-2,
                       tol = 1e-4, zeta = NULL, ...) {
  fit <- fit.additive(y=dat$y, x=dat$x,
                      family="gaussian", method = "tf",
                      lambda.max = lambda.max, lambda.min.ratio = lambda.min.ratio,
                      tol = tol, zeta = zeta,...)

  fit.vals <- apply(fit$f_hat, c(1,3),sum)
  fit.vals <- fit.vals + fit$intercept

  # Evaluate the MSE_true
  mse.true <- colMeans((fit.vals - dat$f0)^2)

  min(mse.true)
}


SimTF <- function(dat, lambda.max = 1, lambda.min.ratio = 1e-2,
                  tol = 1e-4, ...) {
  require(GSAM)

  # lambda.max = 1; lambda.min.ratio = 1e-2
  # tol = 1e-4;max.iter = 300

  # First we fit the lambda, lambda^2 models
  mse0 <- oracle.TF(dat,lambda.max = lambda.max,
                    lambda.min.ratio = lambda.min.ratio,
                    tol = tol, zeta = NULL, ...)

  zeta.len <- 10
  mse.zeta <- numeric(zeta.len)
  zeta.seq <- seq(1e-3, 1-1e-3, length = zeta.len)
  for(i in 1:zeta.len) {
    #print(i)
    mse.zeta[i] <- oracle.TF(dat,lambda.max = lambda.max,
                             lambda.min.ratio = lambda.min.ratio,
                             tol = tol, zeta = zeta.seq[i],...)

  }

  mse1 <- min(mse.zeta)


  return(c("mse0" = mse0,
              "mse1" = mse1))
}


