# This file consists of models we will use in the simulations.
# I.e. we will define the simulation functions here.


# First we define step functions.
func.step <- function(x, knots = c(0), vals = c(-1,1)){
  # knots: The position of the k knots
  # vals: The (k+1) vector of values defining value of step.
  # x: The points at which we wish to evaluate the step function
  myf <- stepfun(knots, vals)
  myf(x)
}

func.lin <- function(x, knots = c(0), vals = c(-1,1)) {
  # knots: The position of the k knots
  # vals: The (k+1) vector of values to be the coefficients of the basis expansion.
  # x: The points at which we wish to evaluate the step function

  bs(x, degree = 1, knots = knots)%*%vals
  # plot(x, bs(x, degree = 1, knots = knots)%*%vals, type = "l")
  # abline(v = c(0), col = "red")
}

# Here we define the 4 functions used in the SPAM paper
# This includes a liner, quadratic, sine and exponential function
func.spam1 <- function(x) {
  -2*sin(2*x)
}

func.spam2 <- function(x) {
  0.8*x^2 - 2.5
}

func.spam3 <- function(x) {
  x - 1/2
}

func.spam4 <- function(x) {
  exp(-0.65*x) - 2.5
}


# Finally we define the hills function
# similar to that from Ryan Tibshirani, trend filtering paper.
func.hills <- function(x, split = 0, vals = c(1, 1, 10),
                       rev = FALSE){

  ans <- rep(NA, length(x))
  if(!rev){
    ans[x<split] <-  vals[1] + sin(vals[2]*x[x<split])

    eps <- (vals[2]/vals[3])*cos(vals[2]*split)/(cos(vals[3]*split))
    delta <- vals[1] + sin(vals[2]*split) - eps*sin(vals[3]*split)

    ans[x>=split] <- delta + eps*sin(vals[3]*x[x>=split])
  } else {
    ans[x>split] <-  vals[1] + sin(vals[2]*x[x>split])

    eps <- (vals[2]/vals[3])*cos(vals[2]*split)/(cos(vals[3]*split))
    delta <- vals[1] + sin(vals[2]*split) - eps*sin(vals[3]*split)

    ans[x<=split] <- delta + eps*sin(vals[3]*x[x<=split])
  }
  ans
}

# Scenario 1: All piecewise constant
scen1 <- function(x) {
  # x: A n*4 matrix for the 4 non-zero functions
  f1 <- function(x){
    func.step(x, knots = c(-2.3, -1.8, -0.5, 1.1),
              vals = c(-3, -2.5, -1, 1, 1.8))
  }
  f2 <- function(x){
    kts <- c(-2, -1, 1, 2)
    vals <- c(3, 1.4, 0, -1.7, -1.8)
    func.step(x, knots = sort(kts),
              vals)
  }
  
  f3 <- function(x){
    func.step(x, knots = c(-1.5, 0.5),
              vals = c(-3.3, 2.5, -1))
  }
  
  f4 <- function(x){
    func.step(x, knots = c(-1.7, -0.4, 1.5, 1.9),
              vals = c(-2.8, 0.3, -1.4, 0.4, 1.8))
  }
  
  cbind(f1(x[,1]), f2(x[,2]), f3(x[,3]), f4(x[,4]))
}

# Scenario 2: All piecewise linear
scen2 <- function(x) {
  # x: A n*4 matrix for the 4 non-zero functions
  f1 <- function(x){
    func.lin(x, knots = c(-2.3, -1.8, -0.5, 1.1),
             vals = c(-3, -2.5, -1, 1, 1.8))
  }
  f2 <- function(x){
    kts <- c(-2, -1, 1, 2)
    vals <- c(3, 1.4, 0, -1.7, -1.8)
    func.lin(x, knots = sort(kts),
             vals)
  }

  f3 <- function(x){
    func.lin(x, knots = c(-1.5, 0.5),
             vals = c(-3.3, 2.5, -1))
  }

  f4 <- function(x){
    func.lin(x, knots = c(-1.7, -0.4, 1.5, 1.9),
             vals = c(-2.8, 0.3, -1.4, 0.4, 1.8))
  }

  cbind(f1(x[,1])+0.5, f2(x[,2])-0.5, f3(x[,3]), f4(x[,4])-1)
}

# Scenario 3: All smooth functions
scen3 <- function(x) {
  # x: A n*4 matrix for the 4 non-zero functions

  cbind(func.spam1(x[,1]), func.spam2(x[,2]), func.spam3(x[,3]), func.spam4(x[,4]))
}

# Scenario 4: Mixture of smooth and
scen4 <- function(x) {
  # x: A n*4 matrix for the 4 non-zero functions
  temp1 <- scen1(x)
  temp2 <- scen2(x)
  cbind(temp1[, 1], temp2[,2], func.spam1(x[,3]), func.spam4(x[,4]))
}


# Scenario 5: Hills functions
scen5 <- function(x) {
  # x: A n*4 matrix for the 4 non-zero functions
  x1 <- func.hills(x[,1], 0, c(1, 1, 12))
  x2 <- func.hills(x[,2], 1, c(1, 2, 8))

  x3 <- func.hills(x[,3], -1, c(0, 3,15), rev = TRUE)

  x4 <- func.hills(x[,4], 1, c(0, 2.5, 10), rev = TRUE)

  cbind(x1,x2,x3,x4)

}



