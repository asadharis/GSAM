#####################################################################################
# # Parkinsons Dataset
# # https://archive.ics.uci.edu/ml/datasets/parkinsons+telemonitoring
# parkinsons <- read.csv("parkinsons.csv")
# # PErform a random suffle of the data.
# set.seed(1)
# parkinsons <- parkinsons[sample(nrow(parkinsons)),]
#
# y <- parkinsons$motor_UPDRS
# x <- parkinsons[, c(2,4,7:22)]
#
# xmat <- apply(x, 2, function(x){(x-min(x))/(max(x) - min(x))})
#
#
#
#
# dat <- list("y" = y,
#             "x" = xmat)
#
# save(dat, file = "parkinsons1.RData")

library(mlbench)
data("BostonHousing")
set.seed(1)
y <- BostonHousing$medv
x <- BostonHousing[, c("crim", "indus", "nox", "rm", "age",
                       "dis", "tax", "ptratio", "b", "lstat")]
xmat <- apply(x, 2, function(x){(x-min(x))/(max(x) - min(x))})
xmat2 <- matrix(runif(ncol(x)*nrow(x)), ncol = ncol(x), nrow = nrow(x))
xmat3 <- apply(xmat,2, sample)

dat <- list("y" = y,
            "x" = cbind(xmat,xmat2, xmat3))

save(dat, file = "BostonHousing.RData")
