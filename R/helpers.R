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
  linear_part <- rowSums(f_hat) + intercept
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
  linear_part <- rowSums(f_hat) + intercept
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
  #n <- length(y)

  #lin_part <- apply(f_hat, 1, sum) + intercept
  lin_part <- rowSums(f_hat) + intercept
  (-y)/(1 + exp(y * lin_part))
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
  #linear_part <- apply(f_hat, 1, sum) + intercept
  linear_part <- rowSums(f_hat) + intercept
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
#   method: "tf" for trendfiltering and "sobolev" for the sobolev norm.
#
# Returns:
#   A list of size 2. The objects are:
#   intercept_new: Updated value of the intercept.
#   ans: Updated value of f_hat, an n*p matrix.
GetZ <- function(f_hat, intercept, step_size,
                 x_mat_ord, ord_mat, k,
                 lambda1, lambda2, y,
                 family = "gaussian",
                 method = "tf", parallel = FALSE,
                 ncores = 8, ...) {

  n <- nrow(f_hat)
  p <- ncol(f_hat)

  if(parallel) {
    #p_grp <- ceiling(min(p/ncores, 100))
    p_grp <- ceiling(p/ncores)
    n_grps <- ceiling(p/p_grp)
    grps <- rep(1:n_grps, each = p_grp)
    grps <- grps[1:p]
  }


  if(family == "gaussian") {
    # Least squares loss
    vector_r <- GetVectorR2(y, f_hat, intercept)
  } else if(family == "binomial") {
    # Logistic loss
    vector_r <- GetVectorR(y, f_hat, intercept)
  }

  intercept_new <- intercept - (step_size * mean(vector_r))

  inter_step <- apply(f_hat, 2, "-", step_size * vector_r)


  ans <- matrix(0, nrow = n, ncol = p)
  if(method == "tf") {
    if(parallel) {
      ans <- tryCatch(temp.func.tf(inter_step, x_mat_ord, ord_mat,
                                   k, lambda1, lambda2, step_size,
                                   n, p, n_grps, p_grp, ...),
                      error = function(e) print(e))


    } else {
      for(i in 1:p) {
        temp_y <- inter_step[ord_mat[, i], i]
        ans[ord_mat[, i], i] <-  solve.prox.tf(temp_y - mean(temp_y),
                                               x_mat_ord[, i],
                                               k = k,
                                               lambda1 = lambda1*step_size,
                                               lambda2 = lambda2*step_size,
                                               ...);
      }
    }
  }

  if(method == "sobolev") {
    if(parallel) {
      ans <- tryCatch(temp.func.spline(inter_step, x_mat_ord, ord_mat,
                                       k, lambda1, lambda2, step_size,
                                       n, p, n_grps, p_grp),
                      error = function(e) print(e))
    } else {
      for(i in 1:p) {
        temp_y <- inter_step[ord_mat[, i], i]
        ans[ord_mat[, i], i] <- cpp_solve_prox(temp_y - mean(temp_y),
                                               x_mat_ord[, i],
                                               lambda1 = lambda1*step_size,
                                               lambda2 = lambda2*step_size,
                                               n = n, n_grid = 1e+3,
                                               lam_tilde_old = 0.5);
      }

    }
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
                       family = "gaussian", method = "tf",
                       parallel = FALSE, ncores = 8,
                       line_search = TRUE,
                       ...) {
  if(!line_search) {
    # Return the next step withou any line search..
    temp_z <- GetZ(f_hat, intercept, step_size,
                   x_mat_ord, ord_mat, k,
                   lambda1, lambda2, y, family = family,
                   method = method, parallel = parallel,
                   ncores = ncores, ...)
    return(temp_z)
  }

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
                   lambda1, lambda2, y, family = family,
                   method = method, parallel = parallel,
                   ncores = ncores, ...)

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


##########################################################
######## Temp functions for parallel computing ###########
##########################################################
temp.func.tf <- function(inter_step, x_mat_ord, ord_mat,
                         k, lambda1, lambda2, step_size,
                         n, p, n_grps, p_grp, ...) {
  foreach(i = 1:n_grps, .combine = cbind) %dopar% {
    ind <- ((i-1)*p_grp + 1):((i-1)*p_grp + p_grp)
    if(i == n_grps) {
      ind <- head(ind, p - ((i-2)*p_grp + p_grp))
    }
    ans_j <- matrix(NA, ncol = length(ind), nrow = n)
    for(j in 1:length(ind)) {
      temp_y <- inter_step[ord_mat[,ind[j]], ind[j]]
      ans_j[ord_mat[, ind[j]],j] <-
        GSAM::solve.prox.tf(temp_y - mean(temp_y),
                            x_mat_ord[, ind[j]],k = k,
                            lambda1 = lambda1*step_size,
                            lambda2 = lambda2*step_size,
                            ...);
    }
    ans_j
  }
}

temp.func.spline <- function(inter_step, x_mat_ord, ord_mat,
                             k, lambda1, lambda2, step_size,
                             n, p,  n_grps, p_grp, ...) {
  foreach(i = 1:n_grps, .combine = cbind) %dopar% {
    ind <- ((i-1)*p_grp + 1):((i-1)*p_grp + p_grp)
    if(i == n_grps) {
      ind <- head(ind, p - ((i-2)*p_grp + p_grp))
    }
    ans_j <- matrix(NA, ncol = length(ind), nrow = n)
    for(j in 1:length(ind)) {
      temp_y <- inter_step[ord_mat[,ind[j]], ind[j]]
      ans_j[ord_mat[, ind[j]],j] <-
        GSAM::cpp_solve_prox(temp_y - mean(temp_y),
                             x_mat_ord[, ind[j]],
                             lambda1 = lambda1*step_size,
                             lambda2 = lambda2*step_size,
                             n = n, n_grid = 1e+3,
                             lam_tilde_old = 0.5);
      # The lam_tilde_old parameter is an experimental feature,
      # does noting at the moment
    }
    ans_j
  }
}
##########################################################
######## Temp functions for parallel computing ###########
##########################################################

