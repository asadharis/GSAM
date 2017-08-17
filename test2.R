library(uniSolve)

set.seed(1)
n <- 500
p <- 3
x <- runif(n, -2.5,2.5)
fx <- ifelse(x>=0, -1, 1)
y <- fx + rnorm(n, sd = 0.1)
xmat <- matrix(x)
ord <- matrix(order(x))

if(p>1 ){
  xmat <- cbind(x, matrix(runif(n*(p-1)), ncol = p-1, nrow = n) )
}
ord <- apply(xmat, 2, order)
y <- rbinom(n, size = 1, prob = 1/(1+exp(-fx)))
#y <- round(1/(1+exp(-fx)))
mp <- mean(y)
y2 <- ifelse(y == 0, -1, 1)

ans <- tf.norm(y, xmat, max.iter = 1000, tol = 1e-20,
                    initpars = NULL, lambda.max = 1, lambda.min.ratio = 1e-4,
                    nlam = 10, k = 0, family = "binomial",
                    initintercept = NULL, step = length(y), alpha = 0.5,
               gamma.par = 0.1)
beep()
plot(xmat[,1], ans$f_hat[,1,6], ylim = c(-2,2))
points(xmat[,1], fx - mean(fx), col = "red", pch = 16, cex = 0.5)

ans2 <- sobolev.norm(y, xmat, family = "binomial",max.iter = 1000, tol = 1e-20,
                     initpars = NULL, initintercept = NULL,
                     lambda.max = 0.02, lambda.min.ratio = 1e-3,
                     nlam = 10, step = length(y), alpha = 0.5,
                     gamma.par = 0.01)

preds <- ans2$f_hat[,1,20] + ans2$intercept[20]

plot(x,1/(1+exp(-fx)))

beep()

plot(xmat[,1], ans2$f_hat[,1,5])
points(xmat[,1], fx - mean(fx), col = "red", pch = 16, cex = 0.5)


plot(c(-2.5,2.5), c(-2, 2), type= "n")
for(i in 70:80){
  points(xmat[,1], ans2$f_hat[,1,i], col = i, pch = 16, cex = 0.5)
}
