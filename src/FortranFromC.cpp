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

// arma::vec l_dot(arma::vec y, arma::vec theta) {
//   // Here we assume that y is in {-1,1}
//
//   return (-1 * y) / (1 + exp(y % theta));
// }

double GetLogistic(arma::vec y,
                   arma::mat f_hat, double intercept) {

  // For the case of logistic loss, we assume that y_i \in {-1,1}
  arma::vec lin_part = sum(f_hat, 1) + intercept;

  return mean(log(1 + exp(-(1*y) % lin_part)));
}

arma::vec GetVectorR(arma::vec y, arma::mat f_hat,
                                 double intercept, int n) {

  // With out notation we calculate the vector r,
  // where r_i is given by
  //  r_i = -(1/n) * l^dot (beta_0 + \sum_j f_ji)
  //
  // For the case of logistic loss, we assume that y_i \in {-1,1}

  arma::vec lin_part = sum(f_hat, 1) + intercept;
  arma::vec temp = (-1 * y)/(1 + exp(y % lin_part));
  return  temp/n;
}

arma::field<mat> GetZ(arma::mat f_hat, double intercept,
                      arma::vec vector_r,
                      int n, int p, double step_size,
                      arma::mat x_mat_ord,
                      arma::umat ord_mat, arma::umat rank_mat,
                      double lambda1, double lambda2) {

  // In notation of our algorithm
  // z is given by prox(x_k - t*nabla(f(x_k)), tg)
  // where the objective is (f + g), t is the step_size and x_k is the k^th
  // iteration.
  //
  // Essentially, this is one step of the proximal gradient descent algorithm for
  // a fixed step size t.

  // First we update the intercept.
  arma::mat intercept_mat(1,1);
  intercept_mat(0,0) = intercept - (step_size * accu(vector_r));

  arma::mat inter_step = f_hat;
  inter_step.each_col() -= step_size * vector_r;

  arma::mat ans(n, p, fill::zeros);
  for(int i = 0; i < p; i++) {
    arma::vec temp_y = inter_step.col(i);

    arma::vec temp_ans = cpp_solve_prox(temp_y.elem(ord_mat.col(i)),
                                        x_mat_ord.col(i), step_size*lambda1/n,
                                        step_size*lambda2/n, n, 500);

    ans.col(i) = temp_ans.elem(rank_mat.col(i));
  }

  arma::field<mat> final_ans(2);
  final_ans(0) = intercept_mat;
  final_ans(1) = ans;

  return final_ans;
}

double Get_f_Z(arma::field<mat> z, arma::vec y) {
  // This function simply calculates f(z),
  // where z is one step of the proximal grad descent as
  // in the previous function and
  // f() is just the logistic loss function.

  // This function is simply a wrapper to make writing the line search
  // algorithm cleaner.
  return GetLogistic(y, z(1), as_scalar(z(0)));
}

double Get_f_hat_Z_X(arma::field<mat> z, arma::field<mat> xk,
                     arma::vec vector_r_xk, double f_xk, double step_size,
                     int p) {
  // Again, only a wrapper to make the line search neater.
  // Args:
  //    z: A field of the one step update of prox grad. Element 1 is intercept and
  //       two is the rest, i.e. f^hats.
  //    xk: The current iteration k, again a field with element 1 k^th iterate of intercept
  //        and element 2 the k^th iterate of the f_js.
  //    vector_r_xk: The vector R, as outputed by function 'GetVectorR' at point xk.
  //    f_xk: The loss function evaluated at point xk
  //    step_size: The value of step size, t.
  //    p: The number of components, used for defining loop.

  // Initialize with f_xk and the terms for the intercept.
  double ans = f_xk +
    accu(vector_r_xk) * as_scalar(z(0) - xk(0)) +
    (1/(2*step_size)) * as_scalar(square(z(0) - xk(0)));

  // Obtain the cross product term.
  double cross_prod_term = accu(vector_r_xk.t() * (z(1) - xk(1)));

  // Obtain the norm term.
  double norm_term = (1/(2*step_size)) * accu(square(z(1) - xk(1)));

  return ans + cross_prod_term + norm_term;
}

arma::field<mat> LineSearch(double alpha, double step_size,
                            arma::vec y,
                            arma::mat f_hat, double intercept,
                            int n, int p,
                            arma::mat x_mat_ord, arma::umat ord_mat,
                            arma::umat rank_mat, double lambda1, double lambda2) {

  arma::vec r_k = GetVectorR(y, f_hat, intercept, n);
  double f_xk = GetLogistic(y, f_hat, intercept);

  // Get all the things together for iterate k, in a field.
  arma::field<mat> xk(2);
  arma::mat intercept_mat(1,1);
  intercept_mat(0,0) = intercept;

  xk(0) = intercept_mat;
  xk(1) = f_hat;

  // Initialize empty things for the algorithm.
  arma::field<mat> temp_z;
  bool convg = false;

  while(!convg) {
    temp_z = GetZ(f_hat, intercept,
                  r_k, n, p, step_size, x_mat_ord,
                  ord_mat, rank_mat,
                  lambda1, lambda2);
    double temp_rhs = Get_f_hat_Z_X(temp_z, xk, r_k, f_xk, step_size, p);
    double temp_lhs =  Get_f_Z(temp_z, y) ;

    if(temp_lhs <= temp_rhs) {
      convg = true;
    } else {
      step_size = alpha * step_size;
    }
  }
  return temp_z;
}



// [[Rcpp::export]]
arma::mat cpp_spp_one(arma::vec y, arma::mat x_ord,
                      arma::umat ord, arma::umat ranks,
                      double lambda1, double lambda2,
                      arma::mat init_fhat, double init_intercept,
                      int n, int p,
                      int max_iter = 100, double tol = 1e-4,
                      double step_size = 1, double alpha = 0.5) {
  // # THis function solves the optimization problem:
  // #
  // # argmin (1/2n) ||y - sum_{j=1}^{p}f_j||_n^2 +
  // #                       \sum_{j=1}^p (lambda1)||f_j||_n +
  // #                                    (lambda2)*sqrt(J(f_j)).
  // #
  // #
  // #
  // # y: The response vector assumed to have elements -1,1.
  // # x_ord: A n*p matrix of ordered covariates.
  // # ord: Matrix of dim n*p giving the orders of elemments for each coavairte.
  // # ranks: Matrix of dim n*p giving the ranks of elements for each coavairte.
  // # lambda1, lambda2: Scalar tuning parameters
  // # init_fhat: Matrix of dim n*p of initial estimates, used for warm starts.
  // # init_intercept: Scalar of initial estimate for intercept.
  // # max_iter: maximum number of iterations for block coordinate descent
  // # tol: Tolerance for algorithm
  // # step: Initial step size choice for the line search algorithm.
  // # alpha: The parameter between 0 and 1 for the line search algorithm.

  int  counter = 0;
  bool converged = false;

  arma::field<mat> old_ans(2);
  arma::mat intercept_mat(1,1);
  intercept_mat(0,0) = init_intercept;
  old_ans(0) = intercept_mat;
  old_ans(1) = init_fhat;

  arma::field<mat> new_ans(2);

  while(counter < max_iter && !converged) {


   new_ans =  LineSearch(alpha, step_size, y, old_ans(1), as_scalar(old_ans(0)),
                         n,p, x_ord, ord, ranks, lambda1, lambda2);


    double temp_res = norm(new_ans(1) - old_ans(1), "fro")/(pow(n*p, 0.5))
    + norm(new_ans(0) - old_ans(0), "fro");
    Rcout << "Criteria: " << temp_res << "\n";
    if(temp_res <= tol) {
      converged = true;
    } else {
      counter = counter + 1;
      old_ans = new_ans;
    }
  }
  return new_ans(1);
}

