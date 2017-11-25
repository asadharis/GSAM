dat <- read.csv("SPAMdat.csv", header = TRUE)

y <- dat[, 58]
x <- dat[,1:57]

ym <- ifelse(y==1, 1, 0)

dat <- list("y" = ym, "x" = x)
save(dat, file = "SPAMdat.RData")
