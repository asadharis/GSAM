# IN this file we have the main functions for the trend-filtering additive models
# We use the c implementation of the package glmgen.

library(glmgen)

# Solves the proximal problem for the sobolev norm
solve.prox.tf <- function(y.ord, x.ord, k = 0, lambda1, lambda2) {
  # We now want to solve the optimization problem.
  # minimize (1/2n) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda1)*||f|| +(lambda2)*TV(f)
  #
  # To do this we first solve the problem
  #   f^hat_lambda2 <- argmin (1/2n) sum(i=1,n) [ y(i) - f(x(i)) ]**2 + (lambda2)*TV(f)

  require(glmgen)

  n <- length(y.ord)
  f_hat <- trendfilter(x = x.ord, y = y.ord,  k = k, family = "gaussian",
              lambda = n*lambda2, thinning = TRUE)$beta[, 1]

  if(length(f_hat) != n){
    # Find how many values are missing
    n.res <- n - length(f_hat)
    # Find where in x.ord we have values too close to each other.
    # Find the smallest 'n.res' values
    ind.s <- order(diff(x.ord))[1:n.res]
    f_hat <- R.utils::insert(f_hat, ind.s, values = f_hat[ind.s])
  }

  #f_hat2 <- genlasso::trendfilter(y = y.ord, pos = x.ord, ord = 0)

  # We return the desired value
  #max(1 - lambda1/sqrt(mean(f_hat^2)), 0)*f_hat
  f_hat

}


#########################################################################
#########################################################################
#########################################################################
############### Begin work on prox grad descent #########################
#########################################################################
#########################################################################
#########################################################################

GetLogistic_tf <- function(y, f_hat, intercept) {

  # For the case of logistic loss, we assume that y_i \in {-1,1}
  lin_part <- apply(f_hat, 1, sum) + intercept

  mean(log(1 + exp(-y * lin_part)))
}

gen_dx1 <- function(n) {
  temp <- diag(rep(-1,n-1))
  temp <- cbind(temp, rep(0,n-1))
  diag(temp[,-1]) <- 1
  temp
}

gen_dx2 <- function(n,x){
  dx1 <- gen_dx1(n)
  dx1[-(n-1),-n] %*% diag(1/diff(sort(x))) %*% dx1
}

gen_dx3 <- function(n,x){
  dx1 <- gen_dx1(n)
  dx1_t <- dx1[-((n-2):(n-1)),-((n-1):n)]
  dx2 <- gen_dx2(n,x)
  dx1_t %*% diag(2/diff(sort(x),2)) %*% dx2
}

gen_d <- function(x,n,k=0){
  if(k==0){
    gen_dx1(n)
  } else if(k==1) {
    gen_dx2(n,x)
  } else if(k==2) {
    gen_dx3(n,x)
  }
}

GetPenalty_tf <- function(f_hat, matD) {
  sum(apply((matD %*% f_hat)^2,2,sum))
}

GetVectorR_tf <- function(y, f_hat, intercept) {
  n <- length(y)

  # With out notation we calculate the vector r,
  # where r_i is given by
  #  r_i = -(1/n) * l^dot (beta_0 + \sum_j f_ji)
  #
  # For the case of logistic loss, we assume that y_i \in {-1,1}

  lin_part <- apply(f_hat, 1, sum) + intercept
  temp <- (-1 * y)/(1 + exp(y * lin_part))
  temp/n
}


GetZ_tf <- function(f_hat, intercept, vector_r, n, p, step_size,
                    x_mat_ord, ord_mat, k,lambda1, lambda2) {

  # In notation of our algorithm
  # z is given by prox(x_k - t*nabla(f(x_k)), tg)
  # where the objective is (f + g), t is the step_size and x_k is the k^th
  # iteration.
  #
  # Essentially, this is one step of the proximal gradient descent algorithm for
  # a fixed step size t.

  # First we update the intercept.

  intercept_new <- intercept - (step_size * sum(vector_r))

  #inter_step <- f_hat
  inter_step <- apply(f_hat, 2, "-", step_size * vector_r)
  ans <- matrix(0, nrow = n, ncol = p)
  for(i in 1:p) {
    temp_y <- inter_step[ord_mat[,i], i]

    ans[ord_mat[, i], i] <-  solve.prox.tf(temp_y, x_mat_ord[,i],
                               k = k, lambda1 = lambda1*step_size/n,
                              lambda2 = lambda2*step_size/n);
  }

  list(intercept_new, ans)
}


Get_f_Z_tf <- function(z, y) {
  # This function simply calculates f(z),
  # where z is one step of the proximal grad descent as
  # in the previous function and
  # f() is just the logistic loss function.

  # This function is simply a wrapper to make writing the line search
  # algorithm cleaner.
  GetLogistic_tf(y, z[[2]], z[[1]])
}

Get_f_hat_Z_X_tf <- function(z, xk, vector_r_xk, f_xk, step_size, p) {
  # Again, only a wrapper to make the line search neater.
  # Args:
  #    z: A field of the one step update of prox grad. Element 1 is intercept and
  #       two is the rest, i.e. f^hats.
  #    xk: The current iteration k, again a field with element 1 k^th iterate of intercept
  #        and element 2 the k^th iterate of the f_js.
  #    vector_r_xk: The vector R, as outputed by function 'GetVectorR' at point xk.
  #    f_xk: The loss function evaluated at point xk
  #    step_size: The value of step size, t.
  #    p: The number of components, used for defining loop.

  # Initialize with f_xk and the terms for the intercept.
  ans <- f_xk +
    sum(vector_r_xk) * (z[[1]] - xk[[1]]) +
    (1/(2*step_size)) * ((z[[1]] - xk[[1]])^2)

  # Obtain the cross product term.
  cross_prod_term <- sum(vector_r_xk %*% (z[[2]] - xk[[2]]));

  # Obtain the norm term.
  norm_term <- (1/(2*step_size)) * sum((z[[2]] - xk[[2]])^2);

  return(ans + cross_prod_term + norm_term)
}

LineSearch_tf <- function(alpha, step_size, y, f_hat, intercept,
                       n, p, x_mat_ord, ord_mat, k, lambda1, lambda2) {

  r_k <-  GetVectorR_tf(y, f_hat, intercept)
  f_xk <- GetLogistic_tf(y, f_hat, intercept);

  # Get all the things together for iterate k, in a field.
  xk <- vector("list", 2)
  xk[[1]] <- intercept
  xk[[2]] <- f_hat

  # Initialize empty things for the algorithm.
  temp_z <- vector("list", 2)
  convg <- FALSE

  # temp_z <- GetZ_tf(f_hat, intercept, r_k, n, p, step_size,
  #                   x_mat_ord, ord_mat, k,lambda1, lambda2)
  # convg <- TRUE
  while(!convg) {
    temp_z <- GetZ_tf(f_hat, intercept, r_k, n, p, step_size,
                      x_mat_ord, ord_mat, k,lambda1, lambda2)
    temp_rhs <- Get_f_hat_Z_X_tf(temp_z, xk, r_k, f_xk, step_size, p);
    temp_lhs <-  Get_f_Z_tf(temp_z, y) ;

    if(temp_lhs <= temp_rhs) {
      convg = TRUE
    } else {
      step_size = alpha * step_size;
    }
  }
  cat("Step size: ", step_size,"\n")
  temp_z
}

#########################################################################
#########################################################################
#########################################################################
################# End work on prox grad descent #########################
#########################################################################
#########################################################################
#########################################################################


tf_logistic_one <- function(y, x_ord, ord, lambda1, lambda2,
                            init_fhat, init_intercept, k=0, n, p,
                            max_iter = 100, tol = 1e-4,
                            step_size = 1, alpha = 0.5) {

  counter <- 1
  converged <- FALSE

  matD <- gen_d(x_ord,n,k=k)
  old_ans <- vector("list", 2);
  old_ans[[1]] <- init_intercept
  old_ans[[2]] <- init_fhat

  new_ans <- vector("list", 2);

  while(counter < max_iter & !converged) {
    new_ans <- LineSearch_tf(alpha, step_size, y, old_ans[[2]], old_ans[[1]],
                             n, p, x_ord, ord, k, lambda1, lambda2)
    # par(ask = TRUE)
    # myte <- old_ans[[2]]
    # plot(x_ord[,1], myte[ord[,1]])

    temp2 <- GetLogistic_tf(y,new_ans[[2]],new_ans[[1]]) + lambda2*GetPenalty_tf(new_ans[[2]], matD)

    cat("Objective function is: ", temp2, "\n")

    temp_res <- sqrt(sum((new_ans[[2]] - old_ans[[2]])^2)/(n*p)) +
      abs(new_ans[[1]] - old_ans[[1]])

    #Rcout << "Criteria: " << temp_res << "\n";
    #cat("criteria is: ", (temp_res), ". Counter: ", counter, "\n")
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

tf_one <- function(y, y.mean, x, x.mean, ord, k = 0, lambda1, lambda2,
                       max.iter = 100, tol = 1e-4, initpars = NULL) {
  # THis function solves the optimization problem:
  #
  # argmin (1/2n) ||y - sum_{j=1}^{p}f_j||_n^2 +
  #                       \sum_{j=1}^p (lambda1)||f_j||_n +
  #                                    (lambda2)*TV(f_j).
  #
  #
  # y: A response vector assumed to be centered
  # y.mean: Mean of uncentered response vector
  # x: A n*p matrix of covariates assumed to be column centered
  # x.mean: Means of uncentered design x.
  # ord: Matrix of dim n*p giving the orders/ranks of each coavairte.
  # k: Order of the trend filtering problem. Possible choices are 0,1,2,3.
  # lambda1, lambda2: Scalar tuning parameters
  # max.iter: maximum number of iterations for block coordinate descent
  # tol: Tolerance for algorithm
  # intpars: Initial parameters, taken as 0 if none provided

  n <- length(y)
  p <- ncol(x)

  if(is.null(initpars)) {
    initpars <- matrix(0, ncol = p, nrow = n)
  }

  old.pars <- initpars

  #p <- 96
  # Begin loop for block coordinate descent
  for(i in 1:max.iter) {
    for(j in 1:p) {
      #print(j)
      res <- y - apply(initpars[, -j], 1, sum)
      initpars[ord[,j], j] <- solve.prox.tf(res[ord[, j]],
                                            x[ord[, j], j], k = k, lambda1, lambda2)
    }

    if(mean((initpars - old.pars)^2) < tol ) {
      return(initpars)
    } else {
      old.pars <- initpars
    }
  }
  return(initpars)
}


tf.norm <- function(y, x, max.iter = 100, tol = 1e-4,
                    initpars = NULL, lambda.max = 1, lambda.min.ratio = 1e-3,
                    nlam = 50, k = 0, family = "binomial",
                    initintercept = NULL, step = 1, alpha = 0.5) {
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


  x <- matrix(x)
  n <- length(y)
  p <- ncol(x)

  if(family == "gaussian") {
    x.mean <- apply(x, 2, mean)
    y.mean <- mean(y)


    x.cen <- scale(x, scale = FALSE)
    y.cen <- y - y.mean
    ord <- apply(x.cen, 2, order)
  } else {
    ord <- apply(x, 2, order)
    ranks <- apply(x, 2, rank)
    x_ord <- apply(x, 2, sort)
  }

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


  if(family == "binomial"){
    y[y==0] <- -1
  }

  ans <- array(0, dim = c(n,p,nlam))
  ans.inters <- numeric(nlam)

  for(i in 1:nlam) {
    print(i)
    if(family == "gaussian") {
      ans[,, i] <- tf_one(y.cen, y.mean, x.cen, x.mean,
                          ord, k= k, lam.seq[i], lam.seq[i]^2,
                          max.iter = max.iter,
                          tol = tol, initpars = initpars)
      initpars <- ans[,, i]
    } else if(family == "binomial") {
      temp <- tf_logistic_one(y, x_ord, ord, lam.seq[i], lam.seq[i]^2,
                              initpars, initintercept, n, p,
                              max_iter = max.iter, tol = tol,
                              step_size = step, alpha = alpha, k=k)
      initpars <- temp$fhat
      initintercept <- temp$intercept

      ans[, , i] <- temp$fhat
      ans.inters[i] <- temp$intercept


    }

  }

  if(family == "gaussian") {
    obj <- list("f_hat" = ans,
                "x.cen" = x.cen,
                "y.cen" = y.cen,
                "x.mean" = x.mean,
                "y.mean" = y.mean,
                "ord" = ord,
                "lam" = lam.seq,
                "k" = k,
                "family" = family)
  } else if(family == "binomial") {
    obj <- list("f_hat" = ans,
                "intercept" = ans.inters,
                "x" = x,
                "ord" = ord,
                "lam" = lam.seq,
                "k" = k,
                "family" = family)
  }

  class(obj) <- "tf"
  return(obj)
}

predict.tf <- function(obj, new.data, type = "function") {

  if(obj$family == "gaussian") {
    new.dat.cen <- scale(new.data, scale = FALSE, center = obj$x.mean)

    nlam <- length(obj$lam)
    p <- dim(obj$f_hat)[2]

    ans <- array(0, dim = c(dim(new.data), nlam) )
    for(i in 1:nlam) {
      for(j in 1:p) {
        ans[,j,i] <- approx(obj$x.cen[,j], obj$f_hat[,j,i], new.dat.cen[, j],
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

  return(ans)
}


