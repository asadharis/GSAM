# Main function for a univariate spline solver
spline_s <- function(y, y.ord, x.ord, theta, lambda){
  # According to the code the default function solves the optimization
  # problem.
  # minimize  (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda)*J(f)
  #
  n <- length(x.ord)
  n.grid <- 500
  .C("callSS_Fortran", y = as.double(y.ord),
     x = as.double(x.ord), sy = as.double(theta),
     lambda = as.double(lambda), n_point = as.integer(n),
     yg = as.double(rep(0,n.grid)), xg = as.double(rep(0,n.grid)), ngrid = as.integer(n.grid),
     iderv = as.integer(2))
}

get.pen <- function(obj) {
  n <- length(obj$xg)
  sum(head(obj$yg, -1)^2*(diff(range(obj$xg))/(n+1)))
}


# Solves the proximal problem for the sobolev norm
solve.prox.spline <- function(y.ord, x.ord, lambda1, lambda2) {
  # We now want to solve the optimization problem.
  # minimize (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda1)*||f|| +(lambda2)*sqrt(J(f))
  #
  # To do this we first solve the prox problem
  #   f^hat_lambda2 <- argmin (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda2)*sqrt(J(f))
  #
  # Which requires us to find the root of thf_hate equation
  # x*sqrt(J)(f^tilde_x) - lambda2/2
  # where
  # f^tilde_x <- argmin (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (x)*J(f)

  require(stats)

   n <- length(y.ord)
  tempf <- function(lam_x) {
    cpp_temp_func(lam_x, y.ord, x.ord,
                  n, n_grid = 1000, lambda2)

    # obj <- cpp_spline(y.ord, x.ord, lam_x, n, 1000)
    # lam_x*sqrt(get.pen(obj)) - lambda2/2
  }
  if(tempf(lambda2*1e+2) < 0) {
    #obj <- spline_s(y, y.ord, x.ord, theta = rep(0, n), lambda2)
    b1 <- cov(x.ord,y.ord)/var(x.ord)
    b0 <- mean(y.ord) - (b1*mean(x.ord))
    fhat <- b0 + (b1*x.ord)
    pen <- 0
  } else {
    lam <- uniroot(tempf, c(lambda2*1e-10,lambda2*1e+2))$root
    f_hat <- cpp_spline(y.ord, x.ord, lam, n, 1000)
    fhat <- f_hat$sy
    pen <- sqrt(get.pen(f_hat))
  }
  temp <- max((1 - lambda1/sqrt(mean(fhat^2))), 0)

  # Finally we return the desired value
  return(list("fhat" = temp*fhat))#,
              #"pen" = temp*pen))
}

# This function calculates the loss function or the logistic loss.
#
# Args:
#   y: The response vector, assumed to take values {-1, 1}
#   f_hat: A n*p matrix where each column is the fitted component function evaluated at the covariate points.
#   intercept: The intercept in of the model.
GetLogistic_spline <- function(y, f_hat, intercept) {
  linear_part <- apply(f_hat, 1, sum) + intercept
  mean(log(1 + exp(-y * linear_part)))
}

GetLS_spline <- function(y, f_hat, intercept) {
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
GetVectorR_spline <- function(y, f_hat, intercept) {
  n <- length(y)

  lin_part <- apply(f_hat, 1, sum) + intercept
  (-y/n)/(1 + exp(y * lin_part))
}
GetVectorR2_spline <- function(y, f_hat, intercept){
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
GetZ_spline <- function(f_hat, intercept, step_size,
                 x_mat_ord, ord_mat,
                 lambda1, lambda2, y,
                 family) {

  n <- nrow(f_hat)
  p <- ncol(f_hat)
  if(family == "gaussian") {
    ############## Continuous #####################
    vector_r <- GetVectorR2_spline(y, f_hat, intercept)
    ############## Continuous #####################

  } else if(family == "binomial") {
    ############## Logistic #####################
    vector_r <- GetVectorR_spline(y, f_hat, intercept)
    ############## Logistic #####################
  }

  intercept_new <- intercept - (step_size * sum(vector_r))

  ans <- matrix(0, nrow = n, ncol = p)
    for(i in 1:p) {
    temp_y <- f_hat[,i] - (step_size*vector_r)

    temp <- solve.prox.spline(temp_y[ord_mat[, i]]- mean(temp_y),
                              x_mat_ord[, i],
                              lambda1 = lambda1*step_size/n,
                              lambda2 = lambda2*step_size/n);
    ans[ord_mat[, i], i] <- temp$fhat
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
LineSearch_spline <- function(alpha, step_size, y, f_hat, intercept,
                       x_mat_ord, ord_mat, lambda1, lambda2,
                       family) {


  if(family == "binomial") {
    ############## Logistic #####################
    loss_val <- GetLogistic_spline(y, f_hat, intercept)
    vector_r <- GetVectorR_spline(y, f_hat, intercept)
    ############## Logistic #####################

  } else if(family == "gaussian") {
    ############## Continuous #####################
    loss_val <- GetLS_spline(y, f_hat, intercept)
    vector_r <- GetVectorR2_spline(y, f_hat, intercept)
    ############## Continuous #####################
  }


  convg <- FALSE
  while(!convg) {
    # Obtain next iteration given step size.
    temp_z <- GetZ_spline(f_hat, intercept, step_size,
                   x_mat_ord, ord_mat,
                   lambda1, lambda2, y, family)

    if(family == "binomial") {
      # Generate the values needed for stopping criteria.
      ############## Logistic #####################
      lhs <- GetLogistic_spline(y, temp_z[[2]], temp_z[[1]])
      ############## Logistic #####################

    } else if(family == "gaussian") {
      ############## Continuous #####################
      lhs <- GetLS_spline(y, temp_z[[2]], temp_z[[1]])
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
  return(c(temp_z, step_size))
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
sobolev_one <- function(y, x_ord, ord, lambda1, lambda2,
                   init_fhat, init_intercept,
                   max_iter = 100, tol = 1e-4,
                   step_size = length(y), alpha = 0.5,
                   family) {

  counter <- 1
  converged <- FALSE
  n <- nrow(x_ord)
  p <- ncol(x_ord)

  # Initialize parameters.
  old_ans <- list(init_intercept, init_fhat)

  while(counter <= max_iter & !converged) {

    new_ans <- LineSearch_spline(alpha, step_size = step_size, y, old_ans[[2]], old_ans[[1]],
                                 x_ord, ord, lambda1, lambda2, family)

    temp_res <- sum((new_ans[[2]] - old_ans[[2]])^2) + (new_ans[[1]] - old_ans[[1]])^2
    temp_res2 <- sum(old_ans[[2]]^2) + (old_ans[[1]])^2
    temp_res <- sqrt(temp_res/temp_res2)

    #cat("Iteration:", counter, ". Res: ", temp_res,"\n")
    if(temp_res <= tol) {
      converged <- TRUE
    } else {
      counter <- counter + 1;
      old_ans <- new_ans;
      step_size <- new_ans[[3]]
    }
  }

  list("fhat" = new_ans[[2]],
       "intercept" = new_ans[[1]],
       "conv" = converged)
}

find.lambdamax <- function(y2, x_mat_ord, ord_mat,
                           gamma, family) {
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

sobolev.norm <- function(y, x, max.iter = 100, tol = 1e-3,
                    initpars = NULL, lambda.max = NULL, lambda.min.ratio = 1e-1,
                    nlam = 30, family,
                    initintercept = NULL, step = length(y), alpha = 0.5,
                    gamma.par = NULL) {
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
  if(is.null(lambda.max)){
    lambda.max <- find.lambdamax(y, x_ord, ord, gamma = gamma.par, family = family)
  }

  lam.seq <- seq(log10(lambda.max), log10(lambda.max*lambda.min.ratio),
                 length = nlam)
  lam.seq <- 10^lam.seq

  ans <- array(0, dim = c(n, p, nlam))
  ans.inters <- numeric(nlam)

  full.ans <- vector("list", length = nlam)

  for(i in 1:nlam) {
    cat("Lambda value: ", i, "\n")
    if(is.null(gamma.par)){
      temp <- sobolev_one(y, x_ord, ord, lam.seq[i], lam.seq[i]^2,
                          init_fhat = initpars, init_intercept = initintercept,
                           max_iter = max.iter, tol = tol,
                          step_size = step, alpha = alpha,
                          family = family)

    } else {
      temp <- sobolev_one(y, x_ord, ord, gamma.par*lam.seq[i],
                          (1-gamma.par)*lam.seq[i],
                          init_fhat = initpars, init_intercept = initintercept,
                           max_iter = max.iter, tol = tol,
                          step_size = step, alpha = alpha,
                          family = family)

    }

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
              "family" = family)

  class(obj) <- "sobolev"
  return(obj)
}





predict.sobolev <- function(obj, new.data, type = "function") {

  # new.data <- matrix(rnorm(n*p), ncol = p)
  if(obj$family == "gaussian") {
    new.dat.cen <- scale(new.data, scale = FALSE, center = obj$x.mean)

    nlam <- length(obj$lam)
    p <- dim(obj$f_hat)[2]

    ans <- array(0, dim = c(dim(new.data), nlam) )
    for(i in 1:nlam) {
      for(j in 1:p) {
        ans[,j,i] <- approx(obj$x.cen[, j], obj$f_hat[,j,i], new.dat.cen[, j],
                            rule = 2)$y
      }
    }
  } else if(obj$family == "binomial") {

    nlam <- length(obj$lam)
    p <- dim(obj$f_hat)[2]

    ans <- array(0, dim = c(dim(new.data), nlam) )
    for(i in 1:nlam) {
      for(j in 1:p) {
        ans[,j,i] <- approx(obj$x[, j], obj$f_hat[,j,i], new.data[, j],
                            rule = 2)$y
      }
    }
  }
  if(type == "function") {
    return(ans)
  } else {
    if(obj$family == "gaussian"){
      return(apply(ans, c(1,3),sum) + obj$y.mean)
    } else {
      temp <- apply(ans, c(1,3),sum)
      temp <- t(apply(temp, 1, "+", obj$intercept))

      return(1/(1 + exp(-temp)))

    }
  }

}

# ssp_one <- function(y, y.mean, x, x.mean, ord, lambda1, lambda2, max.iter = 100, tol = 1e-4,
#                     initpars = NULL) {
#   # THis function solves the optimization problem:
#   #
#   # argmin (1/2n) ||y - sum_{j=1}^{p}f_j||_n^2 +
#   #                       \sum_{j=1}^p (lambda1)||f_j||_n +
#   #                                    (lambda2)*sqrt(J(f_j)).
#   #
#   #
#   #
#   # y: A response vector assumed to be centered
#   # y.mean: Mean of uncentered response vector
#   # x: A n*p matrix of covariates assumed to be column centered
#   # x.mean: Means of uncentered design x.
#   # ord: Matrix of dim n*p giving the orders/ranks of each coavairte.
#   # lambda1, lambda2: Scalar tuning parameters
#   # max.iter: maximum number of iterations for block coordinate descent
#   # tol: Tolerance for algorithm
#   # intpars: Initial parameters, taken as 0 if none provided
#
#   n <- length(y)
#   p <- ncol(x)
#
#   if(is.null(initpars)) {
#     initpars <- matrix(0, ncol = p, nrow = n)
#   }
#
#   old.pars <- initpars
#
#   # Begin loop for block coordinate descent
#   for(i in 1:max.iter) {
#     for(j in 1:p) {
#       res <- y - apply(initpars[, -j],1,sum)
#       initpars[ord[,j], j] <- solve.prox(res, res[ord[,j]], x[ord[,j],j], rep(0, n), lambda1, lambda2)
#     }
#     if(mean((initpars - old.pars)^2) < tol ) {
#       return(initpars)
#     } else {
#       old.pars <- initpars
#     }
#   }
#   return(initpars)
# }
#
#
#
# sobolev.norm <- function(y, x, family = "binomial", max.iter = 100, tol = 1e-4,
#                          initpars = NULL, initintercept = NULL,
#                          lambda.max = 1,
#                          lambda.min.ratio = 1e-3,
#                          nlam = 50, step = length(y), alpha = 0.5,
#                          gamma.par = NULL) {
#   # This function solves the sobolev norm prooblem
#   # In this function we will solve the optmization problem for a
#   # sequence of lambda values using warm starts.
#   # We solve the above lambda values for the lambda1 = lambda and
#   # lambda2 = lambda^2.
#   # set.seed(1)
#   # n <- 200
#   # p <- 4
#   # y = rbinom(n, size = 1, prob = 0.5); x = matrix(rnorm(n*p), ncol = p);
#   # max.iter = 1e+4; tol = 1e-4;
#   # initpars = NULL;
#   # initintercept = NULL;
#   # lambda.max = 1;
#   # lambda.min.ratio = 1e-6;
#   # nlam = 50
#   # step = 1
#   # alpha = 0.5
#   # family = "binomial"
#
#   if(family == "binomial"){
#     y[y==0] <- -1
#   }
#
#   n <- length(y)
#   p <- ncol(x)
#
#   if(family == "gaussian") {
#     x.mean <- apply(x, 2, mean)
#     y.mean <- mean(y)
#
#
#     x.cen <- scale(x, scale = FALSE)
#     y.cen <- y - y.mean
#     ord <- apply(x.cen, 2, order)
#   } else {
#     ord <- apply(x, 2, order) - 1
#     ranks <- apply(x, 2, rank) - 1
#     x_ord <- apply(x, 2, sort)
#   }
#
#   lam.seq <- seq(log10(lambda.max), log10(lambda.max*lambda.min.ratio),
#                  length = nlam)
#   lam.seq <- 10^lam.seq
#
#   if(is.null(initpars)) {
#     initpars <- matrix(0, ncol = p, nrow = n)
#   }
#
#   if(is.null(initintercept) & family == "binomial") {
#     mp <- mean(y)
#     initintercept <- log(mp/(1-mp)) ;
#   }
#
#   if(is.null(initintercept) & family == "gaussian") {
#     initintercept = 0;
#   }
#   ans <- array(0, dim = c(n,p,nlam))
#   ans.inters <- numeric(nlam)
#
#   for(i in 1:nlam) {
#     if(family == "gaussian") {
#       ans[,, i] <- ssp_one(y.cen, y.mean, x.cen, x.mean,
#                            ord, lam.seq[i], lam.seq[i],
#                            max.iter = max.iter,
#                            tol = tol, initpars = initpars)
#       initpars <- ans[,, i]
#     } else if(family == "binomial") {
#
#       if(is.null(gamma.par)){
#         cat("Num lambda: ", i, "\n")
#         temp <- cpp_spp_one(y, x_ord, ord, ranks,
#                             lambda1 = lam.seq[i], lambda2 = lam.seq[i]^2,
#                             initpars, initintercept, n, p,
#                             max_iter = max.iter,  tol = tol,
#                             step_size = step, alpha = alpha)
#       } else {
#         cat("Num lambda: ", i, "\n")
#         temp <- cpp_spp_one(y, x_ord, ord, ranks,
#                             lambda1 = gamma.par*lam.seq[i], lambda2 = (1-gamma.par)*lam.seq[i],
#                             initpars, initintercept, n, p,
#                             max_iter = max.iter,  tol = tol,
#                             step_size = step, alpha = alpha)
#       }
#
#
#       ans[,,i] <- temp$fhat
#       ans.inters[i] <- temp$intercept
#       initpars <- temp$fhat
#       initintercept <- temp$intercept
#
#     }
#
#   }
#
#   if(family == "gaussian") {
#     obj <- list("f_hat" = ans,
#                 "x.cen" = x.cen,
#                 "y.cen" = y.cen,
#                 "x.mean" = x.mean,
#                 "y.mean" = y.mean,
#                 "ord" = ord,
#                 "lam" = lam.seq,
#                 "family" = family)
#   } else if(family == "binomial") {
#     obj <- list("f_hat" = ans,
#                 "intercept" = ans.inters,
#                 "x" = x,
#                 "ord" = ord,
#                 "lam" = lam.seq,
#                 "family" = family)
#   }
#
#   class(obj) <- "sobolev"
#   return(obj)
# }
#


