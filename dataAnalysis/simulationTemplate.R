
simulation.template <- function(x, y,
                                folds, method, ...) {

  p <- ncol(x)
  n <- nrow(x)
  for(i in 1:length(unique(folds))){
    temp.ind <- which(folds == i)
    x.train <- x[temp.ind,]
    y.train <- y[temp.ind]

    temp.x <- scale(x.train)
    xbar <- attributes(temp.x)$'scaled:center'
    x.sd <- attributes(temp.x)$'scaled:scale'

    mod <- method(x = temp.x, y = y.train, family = "gaussian")

    temp.new.x <- scale(x[-temp.ind,], center = xbar, scale = x.sd)

    preds.te <- predict(mod, newx = temp.new.x)
    preds.tr <- predict(mod, newx = temp.x)
    mse.te <- apply((preds.te - y[-temp.ind])^2, 2, mean)
    mse.tr <- apply((preds.tr - y[temp.ind])^2, 2, mean)
    rel.opt <- (mse.te - mse.tr)/mse.tr
  }

}
