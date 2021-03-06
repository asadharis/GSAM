source("Models.R")
# The function to generate different data settings.
#
# Args:
#   n: The sample size.
#   p: The number of covariates
#   noise.var: The variance of the error.
#   scenario: The function used to specify the simulation setting.
#   x.seed: Depecrated, now the seed for x will depend on simulation number.
#
# Returns:
# A list with the following objects:
#   y, y.test, y.val: The data y values
#   f0: The true values of the function at the points x.
#   x, x.test, x.val: The vector of x-values at which we fit our model.
#   FUN: The function equal to the scenario function used.

GenerateData <- function(n = 150, p = 4, noise.var = 1, seed = 1,
                         scenario = scen1, x.seed = 0) {

  #n = 50; p = 6; SNR = 10; seed = 1
  #scenario = scen5
  # Generate a fixed design matrix
  x.seed <- seed %% 50
  #x.seed <- 1
  set.seed(x.seed)
  x <- matrix(runif(n * p, min = -2.5, 2.5), ncol = p, nrow = n)
  y.no.noise <- rowSums(scenario(x[, 1:4]))  # Function values at points x.
  #noise.var <- var(y.no.noise)/SNR  # Variance of error terms.


  # Generate validation and test sets of the same size. Approx half of n
  n.test <- floor(n/2)
  x.test <- matrix(runif(n.test * p, min = -2.5, max = 2.5),
                   ncol = p, nrow = n.test)
  x.val <- matrix(runif(n.test * p, min = -2.5, max = 2.5),
                  ncol = p, nrow = n.test)

  # Now we set the seed to a desired dataset.
  set.seed(seed)
  # Generate response.
  y <- y.no.noise + rnorm(n, sd = sqrt(noise.var))
  y.test <- rowSums(scenario(x.test[, 1:4])) + rnorm(n.test, sd = sqrt(noise.var))
  y.val <- rowSums(scenario(x.val[, 1:4]))+ rnorm(n.test, sd = sqrt(noise.var))

  return(list("x" = x, "y" = y,
              "x.test" = x.test, "y.test" = y.test,
              "x.val" = x.val, "y.val" = y.val,
              "f0" = y.no.noise, "FUN" = scenario))
}
