
# This function fits the SSP penalty and returns a single
# value: the value of the oracle MSE.
oracle.SPLINE <- function(dat,lambda.max = 1, lambda.min.ratio = 1e-2,
                      tol = 1e-4, zeta = NULL, ...) {
  fit <- fit.additive(y=dat$y, x=dat$x,
                      family="gaussian", method = "sobolev",
                      lambda.max = lambda.max, lambda.min.ratio = lambda.min.ratio,
                      tol = tol, zeta = zeta,...)

  fit.vals <- apply(fit$f_hat, c(1,3),sum)
  fit.vals <- fit.vals + fit$intercept

  # Evaluate the MSE_true
  mse.true <- colMeans((fit.vals - dat$f0)^2)

  min(mse.true)
}


SimSPLINE <- function(dat, lambda.max = 1, lambda.min.ratio = 1e-2,
                  tol = 1e-4, ...) {
  require(GSAM)

  # lambda.max = NULL; lambda.min.ratio = 1e-2
  # tol = 1e-4;max.iter = 300

  # First we fit the lambda, lambda^2 models
  mse0 <- oracle.SPLINE(dat,lambda.max = lambda.max,
                    lambda.min.ratio = lambda.min.ratio,
                    tol = tol, zeta = NULL, ...)

  zeta.len <- 10

  zeta.seq <- seq(1e-3, 1-1e-5, length = zeta.len)
  mse.zeta <- tryCatch(ssp_par(dat,lambda.max ,
                                 lambda.min.ratio ,
                                 tol, zeta.seq = zeta.seq,...),
                         error = function(e) print(e))

  mse1 <- min(mse.zeta)


  return(c("mse0" = mse0,
              "mse1" = mse1))
}

##########################################################
######## Temp function for parallel computing ############
##########################################################

ssp_par <- function(dat,lambda.max ,
                    lambda.min.ratio ,
                    tol, zeta.seq = zeta.seq,...) {
  n.zeta <- length(zeta.seq)
  foreach(i = 1:n.zeta, .combine = c) %dopar% {
    fit <- GSAM::fit.additive(y=dat$y, x=dat$x,
                        family="gaussian", method = "sobolev",
                        lambda.max = lambda.max,
                        lambda.min.ratio = lambda.min.ratio,
                        tol = tol, zeta = zeta.seq[i],...)

    fit.vals <- apply(fit$f_hat, c(1,3),sum)
    fit.vals <- fit.vals + fit$intercept

    # Evaluate the MSE_true
    mse.true <- colMeans((fit.vals - dat$f0)^2)

    min(mse.true)
  }
}
