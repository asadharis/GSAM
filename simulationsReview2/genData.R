# A main funtion to generate data for
# timing analysis simulations.

gen_data <- function(seed = 1, n = 300, nsim = 1,
                     sd.noise = 0.1,
                     FUN_f = function(x){0.5*sin(1.5*pi*x)}) {
  set.seed(seed)
  # Generate x matrices for iterations (sorted).
  x_mat <- matrix(runif(n*nsim, min = -1, max = 1),
                  ncol = nsim, nrow = n)
  x_mat <- apply(x_mat, 2, sort)

  # Generate response vectors.
  y_mat <- FUN_f(x_mat) + rnorm(n*nsim, sd = sd.noise)

  return(list("x_mat" = x_mat, "y_mat" = y_mat))
}

