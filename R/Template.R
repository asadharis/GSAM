
# This function fits sparse additive models for one
# value of the tuning parameters lambda1, lambda2. This function
# solves the optimization problem using a proximal gradient descent
# algorithm. The nice feature is that a general form of the
# algorithm can be applied to both least-squares or logistic loss, or
# really any loss function. Properties of the least-squares loss
# allow us to use a block coordinate descent algorithm instead.
# This is a general function which can be used for either
# trend-filtering or the sobolev norm.
#
# Args:
#   y: Response vector of size n.
#   x_ord: The ordered/sorted design matrix of size n*p.
#   ord: A n*p matrix of the ordering for the covariate matrix x.
#   lambda1, lambda2: Tuning parameters.
#   init_fhat: Initial value of the estimated functions.
#   init_intercept: Initial value for the intercept.
#   k: order of the TF.
#   max_iter: Maximum iterations of the algorithm.
#   tol: Tolerance for stopping criteria.
#   step_size: Initial step size for the line search algorithm.
#   alpha: Tuning parameter for the line search algorithm.
#   family: "gaussian" for least squares loss and "binomial" for logistic loss.
#   method: Charecter vector indicating method to use. "tf" for trend-filtering
#           and "sobolev" for the sobolev norm smoothness penalty.
#
# Returns:
#   A list with the following components:
#   fhat: The estimate of the component functions, an n*p matrix.
#   intercept: The estimate of the intercept term, a scalar.
#   conv: A boolean indicator of convergence status.
proxGrad_one <- function(y, x_ord, ord, lambda1, lambda2,
                         init_fhat, init_intercept, k=0,
                         max_iter = 100, tol = 1e-4,
                         step_size = 1, alpha = 0.5,
                         family = "gaussian", method = "tf", ...) {

  # Initialize some objects.
  counter <- 1
  converged <- FALSE
  n <- nrow(x_ord)
  p <- ncol(x_ord)

  # Initialize parameters using warm starts.
  old_ans <- vector("list", 2);
  old_ans[[1]] <- init_intercept
  old_ans[[2]] <- init_fhat

  while(counter < max_iter & !converged) {
    new_ans <- LineSearch(alpha, step_size, y, old_ans[[2]], old_ans[[1]],
                          x_ord, ord, k, lambda1, lambda2, family,
                          method = method, ...)

    temp_res1 <- mean((new_ans[[2]] - old_ans[[2]])^2) + (new_ans[[1]] - old_ans[[1]])^2
    temp_res2 <- mean((old_ans[[2]])^2) + (old_ans[[1]])^2

    # Compare the relative change in parameter values and compare to
    # tolerance to check convergence. Else update iteration.
    if(temp_res1/(temp_res2+1e-30) <= tol) {
      converged <- TRUE
    } else {
      counter <- counter + 1;
      old_ans <- new_ans;
    }
  }

  list("fhat" = new_ans[[2]],
       "intercept" = new_ans[[1]],
       "conv" = converged)
}


# This function performs sparse additive trend-filtering for one
# value of the tuning parameters for the case of a least squares loss.
# This is done via a block-coordinate descent method. This function is
# used when family="gaussian".
#
# Args:
#   y: Response vector of size n.
#   x: The design matrix of size n*p.
#   ord: A n*p matrix of the ordering for the covariate matrix x.
#   lambda1, lambda2: Tuning parameters.
#   init_fhat: Initial value of the estimated functions.
#   k: order of the Trend filtering problem, ignored for sobolev.
#   max_iter: Maximum iterations of the algorithm.
#   tol: Tolerance for stopping criteria.
#   method: Charecter vector indicating method to use. "tf" for trend-filtering
#           and "sobolev" for the sobolev norm smoothness penalty.
#
# Returns:
#   A list with the following components:
#   fhat: The estimate of the component functions, an n*p matrix.
#   intercept: The estimate of the intercept term, a scalar.
#   conv: A boolean indicator of convergence status.

blockCoord_one <- function(y, x, ord,init_fhat, k=0,
                           max_iter = 100, tol = 1e-4,
                           lambda1, lambda2, method = "tf",...) {

  n <- length(y)
  p <- ncol(x)

  counter <- 1
  converged <- FALSE

  old.fhat <- init_fhat

  # Begin loop for block coordinate descent
  while(counter < max_iter & !converged) {
    for(j in 1:p) {
      res <- y - apply(init_fhat[, -j], 1, sum)- mean(y)
      if(method == "tf") {
        init_fhat[ord[,j], j] <- solve.prox.tf(res[ord[, j]] - mean(res),
                                               x[ord[, j], j],
                                               k = k, lambda1, lambda2, ...)
      } else if(method == "sobolev"){
        init_fhat[ord[,j], j] <- solve.prox.spline(res[ord[,j]] - mean(res),
                                                   x[ord[, j], j], lambda1, lambda2)
      }
    }

    temp_res1 <- mean((init_fhat - old.fhat)^2)
    temp_res2 <- mean((old.fhat)^2)

    #cat("Iteration: ", counter,", Tol:", temp_res1/(temp_res2+1e-30),"\n")
    # Compare the relative change in parameter values and compare to
    # tolerance to check convergence. Else update iteration.

    if(temp_res1/(temp_res2+1e-30) <= tol) {
      converged <- TRUE
    } else {
      old.fhat <- init_fhat
      counter <- counter+1
    }
  }
  list("fhat" = init_fhat,
       "intercept" = mean(y),
       "conv" = converged)
}




# This function solves the full problem for a sequence of lambda
# values.
#
# Args:
#   y: Response vector of size n.
#   x: An n*p design matrix.
#   max.iter: Maximum number of iterations for the algorithm
#   tol: The tolerance for the algorithm
#   initpars: Initial parameter values, defaults to NULL f^hat = 0.
#   lambda.max: maximum lambda value
#   lambda.min.ratio: Ratio between largest and smallest lambda value.
#   nlam: Number of lambda values.
#   zeta: If NULL (default) the lambda_1 = lambda and lambda2 = lambda^2. Otherwise
#         a double in [0,1] so that lambda_1 = zeta*lambda and lambda_2 = (1-zeta)*lambda.
#   k: Order of trend-filter
#   family: "Gaussian" for least squares norm, and "binomial" for logistic loss.
#   ininitintercept: Initial value for the intercept term.
#   step: Step size for the proximal graident descent algorithm.
#   alpha: Alpha value for the line search algorithm used in proximal gradient descent.
#   coord.desc: A boolean, indicating if block coordinate descent
#               should be used for family=="gaussian"
#   method: Charecter vector indicating method to use. "tf" for trend-filtering
#           and "sobolev" for the sobolev norm smoothness penalty.

fit.additive <- function(y, x, max.iter = 100, tol = 1e-4,
                    initpars = NULL, lambda.max = 1, lambda.min.ratio = 1e-3,
                    nlam = 50, zeta = NULL,
                    k = 0, family = "binomial",
                    initintercept = NULL, step = length(y), alpha = 0.5,
                    coord.desc = TRUE, method = "tf", ...) {

  # Initialize some values.
  x <- as.matrix(x)
  n <- length(y)
  p <- ncol(x)

  # Generate matrix of orders and ranks.
  ord <- apply(x, 2, order)
  ranks <- apply(x, 2, rank)
  x_ord <- apply(x, 2, sort)

  lam.seq <- seq(log10(lambda.max), log10(lambda.max*lambda.min.ratio),
                 length = nlam)
  lam.seq <- 10^lam.seq

  if(is.null(zeta)) {
    lam1.seq <- lam.seq
    lam2.seq <- lam.seq^2
  } else {
    if(zeta > 1 | zeta < 0) {
      stop("zeta must be within [0,1]")
    }
    lam1.seq <- lam.seq * zeta
    lam2.seq <- lam.seq * (1 - zeta)
  }

  if(is.null(initpars)) {
    initpars <- matrix(0, ncol = p, nrow = n)
  }

  if(is.null(initintercept) & family == "binomial") {
    mp <- mean(y)
    initintercept <- log(mp/(1-mp))
  }
  if(is.null(initintercept) & family == "gaussian") {
    initintercept <- mean(y)
  }

  if(family == "binomial"){
    y[y==0] <- -1
  }

  ans <- array(0, dim = c(n, p, nlam))
  ans.inters <- numeric(nlam)
  conv.vec <- c()

  for(i in 1:nlam) {
    #cat("Lambda value: ", i, "\n")
    if(family=="gaussian" & coord.desc) {
      temp <- blockCoord_one(y, x, ord,init_fhat = initpars,
                             k=k, max_iter = max.iter, tol = tol,
                             lambda1 = lam1.seq[i],lambda2 = lam2.seq[i],
                             method = method,...)
    } else {
      #print(i)
      temp <- proxGrad_one(y, x_ord, ord, lam1.seq[i], lam2.seq[i],
                     init_fhat = initpars, init_intercept = initintercept,
                     k=k, max_iter = max.iter, tol = tol,
                     step_size = step, alpha = alpha,
                     family = family, method = method,...)

    }

    ans[, , i] <- temp$fhat
    ans.inters[i] <- temp$intercept

    initintercept <- ans.inters[i]
    initpars <- ans[, , i]
    conv.vec <- c(conv.vec, temp$conv)
  }

  obj <- list("f_hat" = ans,
              "intercept" = ans.inters,
              "x" = x,
              "ord" = ord,
              "lam" = lam.seq,
              "k" = k,
              "family" = family,
              "conv" = conv.vec)

  class(obj) <- "add_mod"
  return(obj)
}


# The Generic predict function our framework based on linear
# interpolation.
#
# Args:
#   obj: A fitted additive model, object of type "add_mod"
#   new.data: A new x_matrix which we wish to fit.
#   type: "function" will return the fitted components f_1,...,f_p.

predict.add_mod <- function(obj, new.data, type = "function") {
  nlam <- length(obj$lam)
  p <- dim(obj$f_hat)[2]

  ans <- array(0, dim = c(dim(new.data), nlam) )
  for(i in 1:nlam) {
    for(j in 1:p) {
      ans[,j,i] <- approx(obj$x[, j], obj$f_hat[,j,i], new.data[, j],
                          rule = 2)$y
    }
  }

  if(type == "function") {
    return(ans)
  } else {
    if(obj$family == "gaussian"){
      temp <- apply(ans, c(1,3),sum)
      temp <- t(apply(temp, 1, "+", obj$intercept))

      return(temp)
    } else {
      temp <- apply(ans, c(1,3),sum)
      temp <- t(apply(temp, 1, "+", obj$intercept))

      return(1/(1 + exp(-temp)))
    }
  }
  return(ans)
}


# Another generic function to find the lambda_max, based on the proximal
# gradient descent algorithm.
find.lambdamax <- function(y2, x_mat_ord, ord_mat,
                          family) {
  n <- length(y2)
  y.bar <- sum(y2==1)/n
  inter <- log(y.bar/(1-y.bar))


  vecR <- (-y2/n)/(1 + exp(y2 * inter))
  lam.max <- n*sqrt(mean(vecR^2))

  f_hat <- matrix(0, ncol = ncol(x_mat_ord), nrow = nrow(x_mat_ord))

  convg <- FALSE
  while(!convg){
    lam <- lam.max*0.9
    if(is.null(gamma)){
      f_hat_new <- GetZ_spline(f_hat, inter, step_size = 1,
                               x_mat_ord, ord_mat,
                               lam, lam^2, y2, family)
    } else {
      f_hat_new <- GetZ_spline(f_hat, inter, step_size = 1,
                               x_mat_ord, ord_mat,
                               gamma*lam, (1-gamma)*lam, y2, family)
    }
    temp_res <- sum((f_hat_new[[2]] - f_hat)^2)
    if(temp_res < 1e-30){
      lam.max <- lam
    } else {
      convg <- TRUE
    }

  }
  return(lam.max)
}




