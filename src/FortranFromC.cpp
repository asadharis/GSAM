// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include<stdlib.h>
#include<iostream>
#include<math.h>
#include<vector>

using namespace arma;
using namespace Rcpp;
using namespace std;

extern "C"{
  void css_(double *h,
	    int* npoint,
	    double *x,
	    double *y,
	    double *wght,
	    double *sy,
	    double *trace,
	    double *diag,
	    double *cv,
	    int *ngrid,
	    double *xg,
	    double *yg,
	    int *job,
	    int *ideriv,
	    int *ierr);

void callSS_Fortran(double *y, double *x, double *sy, double *lambda, int *n_point,
                    double *yg, double *xg, int *ngrid,
                    int *iderv)
{
  double lam = (*lambda);
  double n = (*n_point);
  double h = log(2*lam); //log lambda
  //double h = log(lam); //log lambda
  int npoint = n; //number of points

  double wght[npoint];
  double trace = 0;
  double diag[npoint];

  for(int i = 0; i < npoint; i++){
    sy[i] = 0;
    wght[i] = 1;
    diag[i] = 0;
  }
  diag[0] = 1;

  double cv = 0;
  //int ngrid = 0;
  //double xg = 0;
  //double yg = 0;

  int job[3]; //Assumes Xs are sorted

  int ideriv = *iderv;
  if(ideriv == 2) {
    job[0] = 0; job[1] = 2; job[2] = 1;
  } else {
    job[0] = 0; job[1] = 0; job[2] = 1;
  }

  int ierr = 0;

  css_(&h,
       &npoint,
       x,
       y,
       wght,
       sy,
       &trace,
       diag,
       &cv,
       ngrid,
       xg,
       yg,
       job,
       &ideriv,
       &ierr);
}
}



// [[Rcpp::export]]
List cpp_spline(arma::vec y_ord, arma::vec  x_ord,
                  double lambda, int n, int n_grid) {

  vec theta(n, fill::zeros);
  vec yg(n_grid, fill::zeros);
  vec xg(n_grid, fill::zeros);

  int iderv = 2;

  callSS_Fortran(y_ord.memptr(), x_ord.memptr(),
                 theta.memptr(), &lambda, &n,
                 yg.memptr(), xg.memptr(), &n_grid, &iderv);

  return List::create(Named("y") = y_ord,
                      Named("x") = x_ord,
                      Named("sy") = theta,
                      Named("lambda") = lambda,
                      Named("n_point") = n,
                      Named("yg") = yg,
                      Named("xg") = xg,
                      Named("ngrid") = n_grid);
}


// [[Rcpp::export]]
arma::vec cpp_spline_raw(arma::vec y_ord, arma::vec  x_ord,
                double lambda, int n, int n_grid) {

  vec theta(n, fill::zeros);
  vec yg(n_grid, fill::zeros);
  vec xg(n_grid, fill::zeros);

  int iderv = 0;

  callSS_Fortran(y_ord.memptr(), x_ord.memptr(),
                 theta.memptr(), &lambda, &n,
                 yg.memptr(), xg.memptr(), &n_grid, &iderv);

  return theta;
}


// [[Rcpp::export]]
double cpp_temp_func(double lambda, arma::vec y_ord, arma::vec x_ord,
                     int n, int n_grid, double lambda2) {

  vec theta(n, fill::zeros);
  vec yg(n_grid, fill::zeros);
  vec xg(n_grid, fill::zeros);

  int iderv = 2;

  callSS_Fortran(y_ord.memptr(), x_ord.memptr(),
                 theta.memptr(), &lambda, &n,
                 yg.memptr(), xg.memptr(), &n_grid, &iderv);

  double diff = xg.max() - xg.min();
  mat ans(1,1);
  ans(0,0) = accu(square(yg.head(n_grid - 1)) * diff/(n_grid + 1));

  return lambda * as_scalar(sqrt(ans)) - lambda2/2.0;
}

// [[Rcpp::export]]
double cpp_uniroot(double lambda_min, double lambda_max,
                   arma::vec y_ord, arma::vec x_ord,
                   int n, int n_grid, double lambda2,
                   double tol) {

  double x_max = lambda_max;
  double x_min = lambda_min;

  if(x_max - x_min <= tol) {
    return (x_max + x_min)/2;
  } else {
    double mid = (x_max + x_min)/2;
    double temp_min = cpp_temp_func(x_min, y_ord, x_ord,
                                    n, n_grid, lambda2);

    double temp = cpp_temp_func(mid, y_ord, x_ord,
                                n, n_grid, lambda2);
    if(temp_min <= 0 && temp >= 0) {
      return cpp_uniroot(x_min, mid, y_ord, x_ord,
                         n, n_grid, lambda2,tol);
    } else if (temp_min <= 0 && temp <= 0) {
      return cpp_uniroot(mid, x_max, y_ord, x_ord,
                         n, n_grid, lambda2, tol);
    } else if (temp_min >= 0 && temp >= 0) {
      return cpp_uniroot(mid, x_max, y_ord, x_ord,
                         n, n_grid, lambda2, tol);
    } else {
      return cpp_uniroot(x_min, mid, y_ord, x_ord,
                         n, n_grid, lambda2, tol);
    }
  }
  return 0;
}



// [[Rcpp::export]]
arma::vec cpp_solve_prox(arma::vec y_ord, arma::vec x_ord,
                         double lambda1, double lambda2,
                         int n, int n_grid) {

  // We now want to solve the optimization problem.
  // minimize (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda1)*||f|| +(lambda2)*sqrt(J(f))
  //
  //   To do this we first solve the prox problem
  //     f^hat_lambda2 <- argmin (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (lambda2)*sqrt(J(f))
  //
  //     Which requires us to find the root of thf_hate equation
  //     x*sqrt(J)(f^tilde_x) - lambda2/2
  //   where
  //     f^tilde_x <- argmin (1/2n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + (x)*J(f)
  //

  arma::vec f_hat;
  double temp_initial = cpp_temp_func(lambda2*100, y_ord, x_ord,
                                      n, n_grid, lambda2);
  if(temp_initial < 0) {
    double b1 = as_scalar(cov(x_ord, y_ord))/as_scalar(var(x_ord));
    double b0 = as_scalar(mean(y_ord) - (b1*mean(x_ord)));
    f_hat = b0 + b1 * x_ord;
  } else {
    double temp_lam = cpp_uniroot(lambda2*1e-10, lambda2*1e+2,
                                  y_ord, x_ord, n, n_grid,
                                  lambda2, 1e-7);
    f_hat = cpp_spline_raw(y_ord, x_ord, temp_lam, n, n_grid);
  }

  double temp_inner = as_scalar(1 - lambda1/mean(square(f_hat)));
  if(temp_inner < 0) {
    return 0*f_hat;
  } else {
    return temp_inner * f_hat;
  }

  return f_hat;

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////// Begin work on proximal gradient descent algorithm.///////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
arma::mat cpp_spp_one(arma::vec y, double y_mean, arma::vec x, double x_mean,
                      arma::mat ord, double lambda1, double lambda2,
                      arma::mat initpars, int n, int p,
                      int max_iter = 100, double tol = 1e-4) {
  // # THis function solves the optimization problem:
  // #
  // # argmin (1/2n) ||y - sum_{j=1}^{p}f_j||_n^2 +
  // #                       \sum_{j=1}^p (lambda1)||f_j||_n +
  // #                                    (lambda2)*sqrt(J(f_j)).
  // #
  // #
  // #
  // # y: A response vector assumed to be centered
  // # y.mean: Mean of uncentered response vector
  // # x: A n*p matrix of covariates assumed to be column centered
  // # x.mean: Means of uncentered design x.
  // # ord: Matrix of dim n*p giving the orders/ranks of each coavairte.
  // # lambda1, lambda2: Scalar tuning parameters
  // # max.iter: maximum number of iterations for block coordinate descent
  // # tol: Tolerance for algorithm
  // # intpars: Initial parameters, taken as 0 if none provided



}

