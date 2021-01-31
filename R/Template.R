
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
#   parallel: A boolean, indicating if proximal gradient descent should
#             be parallelized. Ignored if coord.desc == TRUE.
#   ncores: If parallel == TRUE, specifies size of cluster.
#   line_search: A boolean, indicator if line search should be used or
#                 fixed step size.
#   FISTA: A boolean, should the accelerated proximal gradient descent be used?
#
#
# Returns:
#   A list with the following components:
#   fhat: The estimate of the component functions, an n*p matrix.
#   intercept: The estimate of the intercept term, a scalar.
#   conv: A boolean indicator of convergence status.
proxGrad_one <- function(y, x_ord, ord, lambda1, lambda2,
                         init_fhat, init_intercept, k=0,
                         max_iter = 100, tol = 1e-4,
                         step_size = 5, alpha = 0.8,
                         family = "gaussian", method = "tf",
                         parallel = FALSE, ncores = 8,
                         line_search = TRUE, FISTA = FALSE,
                         ...) {

  # Initialize some objects.
  counter <- 1
  converged <- FALSE
  n <- nrow(x_ord)
  p <- ncol(x_ord)

  # Initialize parameters using warm starts.
  old_ans <- vector("list", 2)
  old_ans[[1]] <- init_intercept
  old_ans[[2]] <- init_fhat

  if(FISTA) {
    # While old_ans will be used to store
    # iterate k, we have another object to store
    # iterate k-1.
    # This is for Nesterov style acceleration.
    old_old_ans <- vector("list", 2)
    old_old_ans[[1]] <- init_intercept
    old_old_ans[[2]] <- init_fhat
  }

  if(line_search == FALSE) {
    # If not doing a line search we need to specify a step size in terms of the
    # Lipchitz constant.
    # Specifcally we use step size ((pa+1)*L)^-1 where L is
    # the lipschitz constant, and pa is the number of non-zero functions on the
    # current iterate.

    if(family == "binomial") {
      L.const <- 1/4
    } else if(family == "gaussian") {
      L.const <- 1
    }
  }

  while(counter < max_iter & !converged) {
    if(line_search == FALSE) {
      sparse_f <- 1 * (colMeans(abs(old_ans[[2]])) != 0)
      step_size <- (L.const * (1 + sum(sparse_f)))^(-1)
    }

    if(FISTA) {
      # Evaluate extrapolation step of fista.
      wk <- counter/(counter+3)
      fista_fhat <- old_ans[[2]] + wk*(old_ans[[2]] - old_old_ans[[2]])
      fista_inter <- old_ans[[1]] + wk*(old_ans[[1]] - old_old_ans[[1]])
      new_ans <- LineSearch(alpha, step_size, y, fista_fhat,
                            fista_inter,
                            x_ord, ord, k, lambda1, lambda2, family,
                            method = method, parallel = parallel,
                            ncores = ncores,
                            line_search = line_search, ...)

    } else {
      new_ans <- LineSearch(alpha, step_size, y,
                            old_ans[[2]], old_ans[[1]],
                            x_ord, ord, k, lambda1, lambda2, family,
                            method = method,parallel = parallel,
                            ncores = ncores,
                            line_search = line_search, ...)
    }


    temp_res1 <- mean((new_ans[[2]] - old_ans[[2]])^2) +
      (new_ans[[1]] - old_ans[[1]])^2
    temp_res2 <- mean((old_ans[[2]])^2) + (old_ans[[1]])^2

    # Compare the relative change in parameter values and compare to
    # tolerance to check convergence. Else update iteration.
    if(temp_res1/(temp_res2+1e-30) <= tol) {
      converged <- TRUE
    } else {
      counter <- counter + 1;
      if(FISTA){
        old_old_ans <- old_ans;
      }
      old_ans <- new_ans;
    }
  }
  if(converged == FALSE) {
    expr <- paste0("Algorithm did not converge for Lambda1 = ", signif(lambda1, 2),
                   " and Lambda2 = ", signif(lambda2,2));
    warning(expr)
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
        init_fhat[ord[,j], j] <- cpp_solve_prox(res[ord[,j]] - mean(res),
                                                   x[ord[, j], j], lambda1, lambda2,
                                                n = n, n_grid = 1e+3,
                                                lam_tilde_old = 0.5)
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
#   family: "gaussian" for least squares norm, and "binomial" for logistic loss.
#   ininitintercept: Initial value for the intercept term.
#   step: Step size for the proximal graident descent algorithm.
#   alpha: Alpha value for the line search algorithm used in proximal gradient descent.
#   coord.desc: A boolean, indicating if block coordinate descent
#               should be used for family=="gaussian"
#   method: Charecter vector indicating method to use. "tf" for trend-filtering
#           and "sobolev" for the sobolev norm smoothness penalty.
#   parallel: A boolean, indicating if proximal gradient descent should
#             be parallelized. Ignored if coord.desc == TRUE.
#   ncores: If parallel == TRUE, specifies size of cluster.
#   line_search: A boolean, indicator if line search should be used or
#                 fixed step size.
#   FISTA: A boolean, should the accelerated proximal gradient descent be used?
#   return_x: A boolean, should the design matrix be returned.
#             for large design matrices we recommend return_x == FALSE.

fit.additive <- function(y, x, max.iter = 100, tol = 1e-4,
                    initpars = NULL, lambda.max = NULL,
                    lambda.min.ratio = 1e-3,
                    nlam = 50, zeta = NULL,
                    k = 0, family = "binomial",
                    initintercept = NULL, step = 1, alpha = 0.8,
                    coord.desc = TRUE, method = "tf",
                    parallel = FALSE, ncores = 8,
                    line_search = FALSE,
                    FISTA = TRUE, verbose = FALSE,
                    return_x = TRUE, ...) {


  method <- tolower(method)
  if(family != "binomial" & family != "gaussian") {
    stop("Currenlty, only 'gaussian' and 'binomial' are acceptable values for
    function parameter: family.")
  }

  if(method != "tf" & method != "sobolev") {
    stop("Currenlty, only 'tf' and 'sobolev' are acceptable values for
    function parameter: method.")
  }
  # Initialize some values.
  x <- as.matrix(x)
  n <- length(y)
  p <- ncol(x)

  # Generate matrix of orders.
  ord <- apply(x, 2, order)
  x_ord <- sapply(1:p, function(i){
    x[ord[,i], i]
  })

  ### Deprecated, no need to calculate ranks,
  ### especially when working with high-dimensional data.
  #ranks <- apply(x, 2, rank)


  if(is.null(initpars)) {
    initpars <- matrix(0, ncol = p, nrow = n)
  }

  if(family == "binomial"){
    y <- factor(y)
    if(length(levels(y))!=2) {
      stop("For binomial family response must be binary or a factor with only two
           levels.")
    }
    levs.y <- levels(y)
    y <- ifelse(y == levels(y)[2], 1, -1)
  }

  if(is.null(initintercept) & family == "binomial") {
    y2 <- y
    y2[y==-1] <- 0
    mp <- mean(y2)
    initintercept <- log(mp/(1-mp))
  }
  if(is.null(initintercept) & family == "gaussian") {
    initintercept <- mean(y)
  }

  # If we are paralleling some computations.
  if(parallel) {
    require(doParallel)
    require(parallel)

    # Begin cluster
    #cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cores = ncores)
    if(verbose) {
      cat("Registering Parallel Backend with", ncores, "Cores.\n")
    }

  }


  if(is.null(lambda.max)) {
    lambda.max <- find.lambdamax(x_ord, ord, k,
    y, family = family,
    method = method, parallel = parallel,
    ncores = ncores, zeta = zeta)
    if(verbose){
      cat("Lambda Max found: ", lambda.max, "\n")
    }

  }


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

  ans <- array(0, dim = c(n, p, nlam))
  ans.inters <- numeric(nlam)
  conv.vec <- c()

  for(i in 1:nlam) {
    if(verbose) {
      cat("Lambda value: ", i, "\n");
    }
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
                     family = family, method = method,
                     parallel = parallel, ncores = ncores,
                     line_search = line_search, FISTA = FISTA,
                     ...)

    }

    ans[, , i] <- temp$fhat
    ans.inters[i] <- temp$intercept

    initintercept <- ans.inters[i]
    initpars <- ans[, , i]
    conv.vec <- c(conv.vec, temp$conv)
  }


  if(return_x == TRUE) {
    obj <- list("f_hat" = ans,
                "intercept" = ans.inters,
                "x" = x,
                "ord" = ord,
                "lam" = lam.seq,
                "k" = k,
                "family" = family,
                "conv" = conv.vec)

  } else {
    obj <- list("f_hat" = ans,
                "intercept" = ans.inters,
                "lam" = lam.seq,
                "k" = k,
                "family" = family,
                "conv" = conv.vec)

  }

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

predict.add_mod <- function(obj, new.data, type = "function",
                            x.mat = obj$x) {
  nlam <- length(obj$lam)
  p <- dim(obj$f_hat)[2]
  obj$x <- x.mat

  ans <- array(0, dim = c(dim(new.data), nlam) )
  for(i in 1:nlam) {
    for(j in 1:p) {
      ans[,j,i] <- suppressWarnings({
        approx(obj$x[, j], obj$f_hat[,j,i], new.data[, j],
               rule = 2)$y
      })
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

plot.add_mod <- function(obj, f_ind = 1, lam_ind = 1,
                         x.mat = obj$x, ...) {
  y <- obj$f_hat[,f_ind,lam_ind]
  plot(x = sort(x.mat[,f_ind]), y=y[order(x.mat[,f_ind])],
       xlab = paste0("X_", f_ind), ylab = paste0("f_", f_ind), ...)
}

# Another generic function to find the lambda_max, based on the proximal
# gradient descent algorithm.
find.lambdamax <- function(x_mat_ord, ord_mat, k,
                           y, family = "gaussian",
                           method = "tf", parallel = FALSE,
                           ncores = 8, zeta = NULL, ...) {
  n <- length(y)
  if(family == "binomial") {
    y.bar <- sum(y==1)/n
    inter <- log(y.bar/(1-y.bar))


    vecR <- (-y)/(1 + exp(y * inter))
    lam.max <- sqrt(mean(vecR^2))

  } else if(family == "gaussian") {
    inter <- mean(y)
    lam.max <- sqrt(mean(y^2))
  }

  f_hat <- matrix(0, ncol = ncol(x_mat_ord), nrow = nrow(x_mat_ord))

  #print(zeta)
  convg <- FALSE
  while(!convg){
    lam <- lam.max*0.9
    if(is.null(zeta)){

      f_hat_new <- GetZ(f_hat, inter, step_size = 1,
                               x_mat_ord, ord_mat, k=k,
                               lam, lam^2, y, family = family,
                        method = method, parallel = parallel,
                        ncores = ncores, ...)
    } else {
      f_hat_new <- GetZ(f_hat, inter, step_size = 1,
                               x_mat_ord, ord_mat,k=k,
                               zeta*lam, (1-zeta)*lam, y,
                        family = family,
                        method = method, parallel = parallel,
                        ncores = ncores, ...)
    }
    temp_res <- mean((f_hat_new[[2]] - f_hat)^2)
    if(temp_res < 1e-30){
      lam.max <- lam
    } else {
      convg <- TRUE
    }

  }
  return(lam.max)
}

