library(GSAM)
library(glmgen)
library(parallel)
library(doParallel)

# We will time test the following parameters
# 1. Parallel vs. Serial
# 2. FISTA vs regular
# 3. Fixed step-size vs line search

# A function to test time.
time.test <- function(n = 100, p = 436, method = "TF" ,
                      parallel = FALSE, FISTA = FALSE,
                      step_fixed = FALSE) {

  # n = 300; p = 8; method = "TF";
  # parallel = FALSE; FISTA = T;
  # step_fixed = FALSE

  set.seed(1)
  x <- matrix(rnorm(n*p), ncol = p, nrow = n)
  y <- sin(2*x[,1]) + x[,2]^2 + rnorm(n, sd = 0.1)

  t.res <- system.time(mod <- fit.additive(y, x, max.iter = 50, tol = 1e-5,
                                           family = "gaussian",lambda.max = 4,
                                           lambda.min.ratio = 1e-3,
                                           coord.desc = FALSE,
                                           method = method, parallel = parallel,
                                          FISTA = FISTA, line_search = !step_fixed,
                                           verbose = FALSE))
  #print(table(mod$conv))
  expr <- paste0("Para = ", parallel, ", FISTA = ", FISTA, ", FixStep = ", step_fixed,
         ": ", round(t.res[3],3), "\n")
  cat(expr)

}

#time.test(parallel = FALSE, FISTA = FALSE, step_fixed = FALSE)

time.test(parallel = FALSE, FISTA = FALSE, step_fixed = TRUE, p= 4000)

#time.test(parallel = FALSE, FISTA = TRUE, step_fixed = FALSE)

#time.test(parallel = FALSE, FISTA = TRUE, step_fixed = TRUE)


#time.test(parallel = TRUE, FISTA = FALSE, step_fixed = FALSE)

time.test(parallel = TRUE, FISTA = FALSE, step_fixed = TRUE, p=4000)

#time.test(parallel = TRUE, FISTA = TRUE, step_fixed = FALSE)

#time.test(parallel = TRUE, FISTA = TRUE, step_fixed = TRUE)
