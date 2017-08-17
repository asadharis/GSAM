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

solve.prox.tf <- function(y.ord, x.ord, k = 0, lambda1, lambda2) {

  # We first solve the problem
  #      f^hat_lambda2 <- argmin (1/2n) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda2)*TV(f)
  # followed by a soft thresholding step.

  require(glmgen)
  n <- length(y.ord)
  #plot(x.ord, y.ord, ylim = c(-1,1))

  # We use n*lambda2 since the default function solves
  #     argmin (1/2) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda2)*TV(f)
  # whereas we would like to have (1/2n) instead of (1/2).
  if(k==2) {
    f_hat <- trendfilter(x = x.ord, y = y.ord,  k = k, family = "gaussian",
                         lambda = n*lambda2, thinning = FALSE,
                         control = trendfilter.control.list(obj_tol = 1e-12,
                                                            max_iter = 600))$beta[, 1]

  } else {
    f_hat <- trendfilter(x = x.ord, y = y.ord,  k = k, family = "gaussian",
                         lambda = n*lambda2, thinning = FALSE)$beta[, 1]

  }

  # We return the desired value after soft-thresolding.
  max(1 - lambda1/sqrt(mean(f_hat^2)), 0)*f_hat
  #f_hat
}


# This function calculates the loss function or the logistic loss.
#
# Args:
#   y: The response vector, assumed to take values {-1, 1}
#   f_hat: A n*p matrix where each column is the fitted component function evaluated at the covariate points.
#   intercept: The intercept in of the model.
GetLogistic <- function(y, f_hat, intercept) {
  linear_part <- apply(f_hat, 1, sum) + intercept
  mean(log(1 + exp(-y * linear_part)))
}

GetLS <- function(y, f_hat, intercept) {
  linear_part <- apply(f_hat, 1, sum) + intercept
  mean((y - linear_part)^2)/2
}

# This function calculates the n-vector 'r', which consists of all derivative information.
# I.e. r_i = -(y_i/n) * {1 + exp[beta + sum(j=1,p) f_j(x_ij) ]}^-1.
# Or more generally we have that
#   r_i = -(1/n) * l^dot (beta_0 + \sum(j=1,p) f_j(x_ij))
#
# Args:
#   y: The response vector y. Assumed to be in {-1,1}
#   f_hat: The n*p matrix of evaluate function values, as in GetLogistic().
#   intercept: The intercept term.
GetVectorR <- function(y, f_hat, intercept) {
  n <- length(y)

  lin_part <- apply(f_hat, 1, sum) + intercept
  (-y/n)/(1 + exp(y * lin_part))
}
GetVectorR2 <- function(y, f_hat, intercept){
  linear_part <- apply(f_hat, 1, sum) + intercept
  n <- length(y)

  -1*(y - linear_part)/n
}
# This function evaluates the next iteration
# for a given set of parameter values and fixed step size.
# In out notation this corresponds to z, where
#   z <- prox_(t g) (theta^k - t \nabla f)
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
GetZ <- function(f_hat, intercept, step_size,
                 x_mat_ord, ord_mat, k,
                 lambda1, lambda2, y,
                 family = "gaussian") {

  n <- nrow(f_hat)
  p <- ncol(f_hat)
  if(family == "gaussian") {
    ############## Continuous #####################
    vector_r <- GetVectorR2(y, f_hat, intercept)
    ############## Continuous #####################

  } else if(family == "binomial") {
    ############## Logistic #####################
    vector_r <- GetVectorR(y, f_hat, intercept)
    ############## Logistic #####################
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
#   For other parameters see function: 'GetZ()'
LineSearch <- function(alpha, step_size, y, f_hat, intercept,
                       x_mat_ord, ord_mat, k, lambda1, lambda2,
                       family = "gaussian") {


  if(family == "binomial") {
    ############## Logistic #####################
    loss_val <- GetLogistic(y, f_hat, intercept)
    vector_r <- GetVectorR(y, f_hat, intercept)
    ############## Logistic #####################

  } else if(family == "gaussian") {
    ############## Continuous #####################
    loss_val <- GetLS(y, f_hat, intercept)
    vector_r <- GetVectorR2(y, f_hat, intercept)
    ############## Continuous #####################

  }


  convg <- FALSE
  while(!convg) {
    # Obtain next iteration given step size.
    temp_z <- GetZ(f_hat, intercept, step_size,
                   x_mat_ord, ord_mat, k,
                   lambda1, lambda2, y, family = family)

    if(family == "binomial") {
      # Generate the values needed for stopping criteria.
      ############## Logistic #####################
      lhs <- GetLogistic(y, temp_z[[2]], temp_z[[1]])
      ############## Logistic #####################

    } else if(family == "gaussian") {
      ############## Continuous #####################
      lhs <- GetLS(y, temp_z[[2]], temp_z[[1]])
      ############## Continuous #####################

    }
    rhs <- loss_val + sum(vector_r)*(temp_z[[1]] - intercept) +
      sum(vector_r %*% (temp_z[[2]] - f_hat)) +
      (1/(2*step_size)) * (sum((temp_z[[2]] - f_hat)^2) + (temp_z[[1]] - intercept)^2)

    # Check for convergence.
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
# value of the tuning parameters lambda1, lambda2.
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
tf_one <- function(y, x_ord, ord, lambda1, lambda2,
                   init_fhat, init_intercept, k=0,
                   max_iter = 100, tol = 1e-4,
                   step_size = 1, alpha = 0.5,
                   family = "gaussian") {

  counter <- 1
  converged <- FALSE
  n <- nrow(x_ord)
  p <- ncol(x_ord)

  # Initialize parameters.
  old_ans <- vector("list", 2);
  old_ans[[1]] <- init_intercept
  old_ans[[2]] <- init_fhat


  while(counter < max_iter & !converged) {
    new_ans <- LineSearch(alpha, step_size, y, old_ans[[2]], old_ans[[1]],
                          x_ord, ord, k, lambda1, lambda2, family)

    temp_res <- sum((new_ans[[2]] - old_ans[[2]])^2)/(n*p) +
      (new_ans[[1]] - old_ans[[1]])^2

    if(temp_res <= tol) {
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


tf.norm <- function(y, x, max.iter = 100, tol = 1e-4,
                    initpars = NULL, lambda.max = 1, lambda.min.ratio = 1e-3,
                    nlam = 50, k = 0, family = "binomial",
                    initintercept = NULL, step = 1, alpha = 0.5,
                    gamma.par = 0.5) {
  # This function solves the trenfiltering/total variation problem.
  # In this function we will solve the optmization problem for a
  # sequence of lambda values using warm starts.
  # We solve the above lambda values for the lambda1 = lambda and
  # lambda2 = lambda^2.

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


  x <- as.matrix(x)
  n <- length(y)
  p <- ncol(x)

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
    temp <- tf_one(y, x_ord, ord, gamma.par*lam.seq[i], (1-gamma.par) * lam.seq[i],
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


