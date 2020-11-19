# In this file we test our functions for binary data.
library(GSAM)
set.seed(1)
n <- 300; p <- 27;
x <- matrix(rnorm(n*p), ncol = p, nrow = n)
y.lin <- 2*sin(2*x[,1]) + (.5*x[,2])^2
y.prob <- 1/(1+exp(-y.lin))
y <- rbinom(n, size = 1, prob = y.prob)
boxplot(y.prob~y)
#y <- 1*(y.prob>=0.5)

mod_tf1 <- fit.additive(y, x, max.iter = 500, tol = 1e-5,
                     family = "binomial", verbose = TRUE,
                     method = "tf", k=1, lambda.max = 0.1,
                     lambda.min.ratio = 1e-4)
mod_spp <- fit.additive(y, x, max.iter = 300, tol = 1e-4,
                     family = "binomial", verbose = TRUE,
                     lambda.max = .1, FISTA = TRUE,
                     lambda.min.ratio = 1e-4,
                     method = "sobolev", k=0, line_search = TRUE)

#plot(mod_tf0, f_ind = 1, lam_ind = 45, type = "l", col = "red")

plot(mod_spp,f_ind = 1, lam_ind = 7, type = "l", col = "red")
points(x[,1], 2*sin(2*x[,1]) - mean(2*sin(2*x[,1])), col= "blue", pch = 16, cex = 0.6)

plot(mod_spp,f_ind = 2, lam_ind = 4, type = "l", col = "red")
points(x[,2], (.5*x[,2])^2 - mean((.5*x[,2])^2), col= "blue", pch = 16, cex = 0.6)

plot(mod_tf1,f_ind = 1, lam_ind = 8, type = "l", col = "red")
#points(x[,2], (.5*x[,2])^2 - mean((.5*x[,2])^2), col= "blue", pch = 16, cex = 0.6)
points(x[,1], 2*sin(2*x[,1]) - mean(2*sin(2*x[,1])), col= "blue", pch = 16, cex = 0.6)
