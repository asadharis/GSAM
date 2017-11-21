library(glmgen)

###############################################################
#
# This file consists of a proximal gradient descent algorithm
# implementation for the trend-filtering problem.
#
###############################################################

# A function for the proximal operator which solves:
#
# minimize (1/2n) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda1)*||f|| +(lambda2)*TV(f)
#
# where TV(f) is the trend filtering penalty.
# Args:
#   y.ord: The response vector, ordered according to the ordering of x.
#   x.ord: The sorted covariate vector.
#   k: The order of the trend-filtering. We consider 0,1,2
#   lambda1,lambda2: The two tuning parameters.
#   ...: Other parameters that can be passed onto the function glmgen::trendfilter
#
# Returns:
#   f_hat: An n-vector, the solution of the optimization problem.
solve.prox.tf <- function(y.ord, x.ord, k = 0, lambda1, lambda2, ...) {

  # We first solve the problem
  #      f^hat_lambda2 <- argmin (1/2n) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda2)*TV(f)
  # followed by a soft thresholding step.

  require(glmgen)
  n <- length(y.ord)

  # We use n*lambda2 since the default function solves
  #     argmin (1/2) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda2)*TV(f)
  # whereas we would like to have (1/2n) instead of (1/2).

  f_hat <- trendfilter(x = x.ord, y = y.ord,  k = k, family = "gaussian",
                       lambda = n*lambda2, thinning = TRUE, ...)$beta[, 1]

  # The ... allow for other parameters to be passed, we normally would like to
  # pass the control object, e.g. the control below increases tolerance and the
  # maximum number of iterations.
  #
  # control = trendfilter.control.list(obj_tol = 1e-12, max_iter = 600)

  # The following procedure is used becaused of the thinning feature
  # of trendfiltering solver. While this is very useful for making the solver
  # more efficient this can lead to resulting vectors being smaller than n.
  # The if{} block below aims to rectify this issue.
  if(length(f_hat) != n){
    # Find how many values are missing
    n.res <- n - length(f_hat)
    # Find where in x.ord we have values too close to each other. And find the smallest 'n.res' values
    ind.s <- order(diff(x.ord))[1:n.res]
    f_hat <- R.utils::insert(f_hat, ind.s, values = f_hat[ind.s])
  }

  # We return the desired value after soft-thresolding.
  max(1 - lambda1/sqrt(mean(f_hat^2)), 0)*f_hat
}


# Evaluates loss function at given values of parameters, for logistic loss.
#
# Args:
#   y: The response vector, assumed to take values {-1, 1}
#   f_hat: A n*p matrix where each column is the fitted component function evaluated at the covariate points.
#   intercept: The intercept in of the model.
#
# Returns:
#   A scalar value of the loss function.
GetLogistic <- function(y, f_hat, intercept) {
  linear_part <- apply(f_hat, 1, sum) + intercept
  mean(log(1 + exp(-y * linear_part)))
}

# Evaluates loss function at given values of parameters, for least squares loss.
#
# Args:
#   Same as function GetLogistic()
#
# Results:
# A scalar value of the loss.
GetLS <- function(y, f_hat, intercept) {
  linear_part <- apply(f_hat, 1, sum) + intercept
  mean((y - linear_part)^2)/2
}

# This function calculates the n-vector 'r', which consists of all derivative information.
# I.e. r_i = -(y_i/n) * {1 + exp[beta + sum(j=1,p) f_j(x_ij) ]}^-1.
# Or more generally we have that
#   r_i = -(1/n) * l^dot (beta_0 + \sum(j=1,p) f_j(x_ij))
# This particular function calulates r for the logistic loss.
#
# Args:
#   y: The response vector y. Assumed to be in {-1,1}
#   f_hat: The n*p matrix of evaluate function values, as in GetLogistic().
#   intercept: The intercept term.
#
# Returns:
#   r: An n-vector, the derivative vector.
GetVectorR <- function(y, f_hat, intercept) {
  n <- length(y)

  lin_part <- apply(f_hat, 1, sum) + intercept
  (-y/n)/(1 + exp(y * lin_part))
}

# Another function for calculating the 'r' vector.
#   r_i = -(1/n) * l^dot (beta_0 + \sum(j=1,p) f_j(x_ij))
# This particular function calulates r for the least squares loss.
#
# Args:
#   y: The response vector y.
#   f_hat: The n*p matrix of evaluate function values, as in GetLogistic() or GetLS().
#   intercept: The intercept term.
#
# Returns:
#   r: An n-vector, the derivative vector.
GetVectorR2 <- function(y, f_hat, intercept){
  linear_part <- apply(f_hat, 1, sum) + intercept
  n <- length(y)

  -1*(y - linear_part)/n
}


# This function evaluates the an iteration
# for the proximal gradient descent method,
# for a given set of parameter values and fixed step size.
# In our notation it calculates z, where
#   z <- prox_(t * g) (theta^k - t \nabla f)
# where theta^k is the current parameter value, t is the step size,
# f and g are the differentiable and non-differentiable parts of the objective
#     f + g
#
# Args:
#   f_hat, intercept: Values of current iterate, theta^k above.
#   step_size: The step size, written as 't' above.
#   x_mat_ord: The n*p desgin matrix with all columns sorted.
#   ord_mat: An n*p matrix for the ordering of the original (unsorted) design matrix.
#   k: The order of the trend-filtering that we will use.
#   lambda1, lambda2: The tuning parameters.
#   y: The response vector to evaluate the derivate.
#   family: Indicating the loss function we are using, "gaussian" for
#           least squares loss and "binomial" for logistc loss.
#
# Returns:
#   A list of size 2. The objects are:
#   intercept_new: Updated value of the intercept.
#   ans: Updated value of f_hat, an n*p matrix.
GetZ <- function(f_hat, intercept, step_size,
                 x_mat_ord, ord_mat, k,
                 lambda1, lambda2, y,
                 family = "gaussian") {

  n <- nrow(f_hat)
  p <- ncol(f_hat)
  if(family == "gaussian") {
    # Least squares loss
    vector_r <- GetVectorR2(y, f_hat, intercept)
  } else if(family == "binomial") {
    # Logistic loss
    vector_r <- GetVectorR(y, f_hat, intercept)
  }

  intercept_new <- intercept - (step_size * sum(vector_r))

  inter_step <- apply(f_hat, 2, "-", step_size * vector_r)

  ans <- matrix(0, nrow = n, ncol = p)
  for(i in 1:p) {
    temp_y <- inter_step[ord_mat[, i], i]

    ans[ord_mat[, i], i] <-  solve.prox.tf(temp_y - mean(temp_y),
                                           x_mat_ord[, i],
                                           k = k, lambda1 = lambda1*step_size/n,
                                           lambda2 = lambda2*step_size/n);
  }
  list(intercept_new, ans)
}

# This function performs the line search algorithm and returns the last
# update of the parameters of the algorithm.
#
# Args:
#   alpha: The parameter for the line search in (0, 1).
#   step_size: An initial value for the step size. Common value is 1.
#   Other parameters: see function, 'GetZ()'
#
# Returns:
#   The next update of the algorithm based on line search step size.
#   The resulting object has the same form as that of the GetZ() function.
LineSearch <- function(alpha, step_size, y, f_hat, intercept,
                       x_mat_ord, ord_mat, k, lambda1, lambda2,
                       family = "gaussian") {

  if(family == "binomial") {
    # Logistic Loss
    loss_val <- GetLogistic(y, f_hat, intercept)
    vector_r <- GetVectorR(y, f_hat, intercept)
  } else if(family == "gaussian") {
    # Least squares loss
    loss_val <- GetLS(y, f_hat, intercept)
    vector_r <- GetVectorR2(y, f_hat, intercept)
  }

  convg <- FALSE
  while(!convg) {
    # Obtain next iteration given step size.
    temp_z <- GetZ(f_hat, intercept, step_size,
                   x_mat_ord, ord_mat, k,
                   lambda1, lambda2, y, family = family)

    # Generate the values needed for stopping criteria.
    if(family == "binomial") {
      lhs <- GetLogistic(y, temp_z[[2]], temp_z[[1]])
    } else if(family == "gaussian") {
      lhs <- GetLS(y, temp_z[[2]], temp_z[[1]])
    }

    rhs <- loss_val + sum(vector_r)*(temp_z[[1]] - intercept) +
      sum(vector_r %*% (temp_z[[2]] - f_hat)) +
      (1/(2*step_size)) * (sum((temp_z[[2]] - f_hat)^2) + (temp_z[[1]] - intercept)^2)

    # Check for convergence else update step size.
    if(lhs <= rhs) {
      convg <- TRUE
    } else {
      step_size <- alpha * step_size
    }
  }

  # Return the next iterate based on the selected step size.
  return(temp_z)
}

# This function performs sparse additive trend-filtering for one
# value of the tuning parameters lambda1, lambda2. This function
# solves the optimization problem using a proximal gradient descent
# algorithm. The nice feature is that a general form of the
# algorithm can be applied to both least-squares or logistic loss, or
# really any loss function. Properties of the least-squares loss
# allow us to use a block coordinate descent algorithm instead.
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
#
# Returns:
#   A list with the following components:
#   fhat: The estimate of the component functions, an n*p matrix.
#   intercept: The estimate of the intercept term, a scalar.
#   conv: A boolean indicator of convergence status.
tf_one <- function(y, x_ord, ord,
                   lambda1, lambda2,
                   init_fhat, init_intercept, k=0,
                   max_iter = 100, tol = 1e-4,
                   step_size = 1, alpha = 0.5,
                   family = "gaussian") {

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
                          x_ord, ord, k, lambda1, lambda2, family)

    temp_res1 <- mean((new_ans[[2]] - old_ans[[2]])^2) + (new_ans[[1]] - old_ans[[1]])^2
    temp_res2 <- mean((old_ans[[2]])^2) + (old_ans[[1]])^2

    # Compare the relative change in parameter values and compare to
    # tolerance to check convergence. Else update iteration.
    if(temp_res1/temp_res2 <= tol) {
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
# This is done via a blocl-coordinate descent method. This function is
# used when family="gaussian".
#
# Args:
#   y: Response vector of size n.
#   x_mat: The design matrix of size n*p.
#   ord: A n*p matrix of the ordering for the covariate matrix x.
#   lambda1, lambda2: Tuning parameters.
#   init_fhat: Initial value of the estimated functions.
#   k: order of the TF.
#   max_iter: Maximum iterations of the algorithm.
#   tol: Tolerance for stopping criteria.
#
# Returns:
#   A list with the following components:
#   fhat: The estimate of the component functions, an n*p matrix.
#   intercept: The estimate of the intercept term, a scalar.
#   conv: A boolean indicator of convergence status.

tf_one_block_coord <- function(y, x_mat, ord,
                               init_fhat, k=0,
                               max_iter = 100, tol = 1e-4,
                               lambda1, lambda2) {

  n <- length(y)
  p <- ncol(x_mat)

  counter <- 1
  converged <- FALSE

  old.fhat <- init_fhat

  # Begin loop for block coordinate descent
  while(counter < max_iter & !converged) {
    for(j in 1:p) {
      res <- y - apply(init_fhat[, -j], 1, sum)
      init_fhat[ord[,j], j] <- solve.prox.tf(res[ord[, j]] - mean(res),
                                             x[ord[, j], j],
                                             k = k, lambda1, lambda2)
    }

    temp_res1 <- mean((init_fhat - old.fhat)^2)
    temp_res2 <- mean((old.fhat)^2)

    # Compare the relative change in parameter values and compare to
    # tolerance to check convergence. Else update iteration.
    if(temp_res1/temp_res2 <= tol) {
      converged <- TRUE
    } else {
      old.fhat <- init_fhat
    }
  }
  list("fhat" = init_fhat,
       "intercept" = mean(y),
       "conv" = converged)
}

# This function solves the trenfiltering/total variation problem.
# In this function we will solve the optmization problem for a
# sequence of lambda values using warm starts.
#
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
#   k: Order of trend-filter
#   family: "Gaussian" for least squares norm, and "binomial" for logistic loss.
#   ininitintercept: Initial value for the intercept term.
#   step: Step size for the proximal graident descent algorithm.
#   alpha: Alpha value for the line search algorithm used in proximal gradient descent.
#   coord.desc: A boolean, indicating if block coordinate descent should be used for family=="gaussian"

tf.norm <- function(y, x, max.iter = 100, tol = 1e-4,
                    initpars = NULL, lambda.max = 1, lambda.min.ratio = 1e-3,
                    nlam = 50, k = 0, family = "binomial",
                    initintercept = NULL, step = 1, alpha = 0.5,
                    coord.desc = TRUE) {


  # n <- 200
  # p <- 100
  # x = matrix(rnorm(n*p), ncol = p);
  # y = sin(3*x[,1]) + rnorm(n,sd = 0.1);
  # max.iter = 100; tol = 1e-4;
  # initpars = NULL;
  # lambda.max = 1;
  # lambda.min.ratio = 1e-3;
  # nlam = 50
  # k <- 0

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

  for(i in 1:nlam) {
    #cat("Lambda value: ", i, "\n")
    if(family=="gaussian" & coord.desc) {
      temp <- tf_one_block_coord(y, x, ord,init_fhat = initpars,
                                 k=k, max_iter = max.iter, tol = tol,
                                 lambda1 = lam.seq[i],lambda2 = lam.seq[i]^2)

      ans[, , i] <- temp$fhat
      ans.inters[i] <- temp$intercept

      initintercept <- ans.inters[i]
      initpars <- ans[, , i]
    }
    temp <- tf_one(y, x_ord, ord, lam.seq[i], lam.seq[i]^2,
                          init_fhat = initpars, init_intercept = initintercept,
                          k=k, max_iter = max.iter, tol = tol,
                          step_size = step, alpha = alpha,
                          family = family)
    ans[, , i] <- temp$fhat
    ans.inters[i] <- temp$intercept

    initintercept <- ans.inters[i]
    initpars <- ans[, , i]
  }

  obj <- list("f_hat" = ans,
              "intercept" = ans.inters,
              "x" = x,
              "ord" = ord,
              "lam" = lam.seq,
              "k" = k,
              "family" = family)

  class(obj) <- "tf"
  return(obj)
}

# The Generic predict function for trend filtering.

predict.tf <- function(obj, new.data, type = "function") {


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


