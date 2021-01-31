# An R file for writing the prox problem for SpAM.

#library(SAM)
library(splines)

# The proximal problem solves the following:
#
# (1/2) * |r - f|_n^2 + lam1|f|_n + lam2P(f)
#
# For spam this problem is equivalent to
# (1/2) * |Qt*r - beta|_n^2 + lam1| beta |_n
#
# where Q is such that X = QR and Qt*Q = I, and
# X is the design matrix for SpAM.
#
# In the case of SpAM, the P(f) is simply a projection onto
# the linear space of basis functions.

# Helper function to generate Q matrix.
Qmat_spam <- function(x, df = 3) {
  x_mat <- ns(x, df = df)
  x_mat_qr <- qr(x_mat)
  x_matQ <- qr.Q(x_mat_qr)
  #x_matR <- qr.R(x_mat_qr)
  return(x_matQ)
}

prox_spam <- function(r, x_matQ, lam1 = 1) {
  r_temp <- crossprod(x_matQ, r)
  r_temp_norm <- sqrt(mean(r_temp^2))

  as.numeric(x_matQ %*% (max((1 - lam1/r_temp_norm),0)*r_temp))
}


