library(uniSolve)
source('Generate_Data.R')
source('Models.R')

seed=1; n = 500;
num.vars = 6; noise.var = 1;
scen.num <- 1
scen = scen1
dat <- GenerateData(seed = seed, n = n, p = num.vars,
                    noise.var = noise.var, scenario = scen)
fit <- fit.additive(y=dat$y, x=dat$x,
                    family="gaussian", method = "sobolev",
                    lambda.max = 0.5, lambda.min.ratio = 1e-4,
                    tol = 5e-5, max.iter = 3000, coord.desc = TRUE)
fit.vals <- apply(fit$f_hat, c(1,3),sum)
fit.vals <- fit.vals + fit$intercept

# Evaluate the MSE_true
mse.true <- colMeans((fit.vals - dat$f0)^2)
plot(mse.true, log = "y")
xout <- seq(-2.5,2.5,length = 500)
xnew <- cbind(xout,xout,xout,xout,xout,xout)
fhat <- predict(fit, xnew)[,,14]
plot(xnew[,4], fhat[,4], type = "l")

#############################

fit2 <- blockCoord_one(dat$y, dat$x, apply(dat$x,2,order),
                       init_fhat = fit$f_hat[,,45], k=0,
                       max_iter = 3000, tol = 1e-5,
                        fit$lam[46], fit$lam[46]^2, method = "sobolev")

fit3 <- proxGrad_one(dat$y, apply(dat$x,2,sort), apply(dat$x,2,order),
                     fit$lam[46], fit$lam[46]^2, init_fhat=fit$f_hat[,,45],
                     init_intercept = fit$intercept[45],
                     k=0,max_iter = 1000, tol = 1e-5,
                     step_size = 500, alpha = 0.5,
                     family = "gaussian", method = "sobolev")
