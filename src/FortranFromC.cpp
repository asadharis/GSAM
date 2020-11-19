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

  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  /////// EDIT TO ORIGINAL SOURCE CODE////////////////////
  /////////////////////////////////////////////////////////
  // Here instead of passing log(lam) we pass
  // log(2*lam). This is because the Fortron routine
  // solves the problem:
  //
  // (1/n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + lambda*J(f)
  //
  // whereas in our paper we have 1/(2n) instead of (1/n).
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
double cpp_temp_func2(double lambda, arma::vec y_ord, arma::vec x_ord,
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

  return as_scalar(sqrt(ans));
}


// [[Rcpp::export]]
double cpp_find_lamdaMax(arma::vec y_ord, arma::vec x_ord,
                       int n, int n_grid, double lam_old,
                       double lambda2) {
  double zero_num = 0.0;
  int iderv = 2;

  vec theta(n, fill::zeros);
  vec yg(n_grid, fill::zeros);
  vec xg(n_grid, fill::zeros);

  callSS_Fortran(y_ord.memptr(), x_ord.memptr(),
                 theta.memptr(), &zero_num, &n,
                 yg.memptr(), xg.memptr(), &n_grid, &iderv);
  arma::mat integ_trap = arma::trapz(xg, square(yg));
  double lam_max = as_scalar(integ_trap);
  double temp;
  temp = cpp_temp_func(lam_max, y_ord, x_ord,
                       n,n_grid,lambda2);
  if(temp < 0) {
    return -1;
  }

  if(lambda2 < lam_max) {
   temp = cpp_temp_func(lambda2, y_ord, x_ord,
                                n,n_grid,lambda2);
    if(temp > 0) {
      //Rcout << "Here in lam2 \n";
      lam_max = lambda2;
    }
  }

  if(lam_old < lam_max) {
  temp = cpp_temp_func(lam_old, y_ord, x_ord,
                                n,n_grid,lambda2);
    if(temp> 0) {
      //Rcout << "Here in lam_old \n";
      lam_max = lam_old;
    }
  }
  return lam_max;
}


double r8_epsilon ( ){
  const double value = 2.220446049250313E-016;

  return value;
}


// [[Rcpp::export]]
double cpp_uniroot3(double a, double b,
                    arma::vec y_ord, arma::vec x_ord,
                    int n, int n_grid, double lambda2,
                    double t){
  double c;
  double d;
  double e;
  double fa;
  double fb;
  double fc;
  double m;
  double macheps;
  double p;
  double q;
  double r;
  double s;
  double sa;
  double sb;
  double tol;
  //
  //  Make local copies of A and B.
  //
  sa = a;
  sb = b;
  fa = cpp_temp_func(sa, y_ord, x_ord, n, n_grid, lambda2);
  fb = cpp_temp_func(sb, y_ord, x_ord, n, n_grid, lambda2);

  c = sa;
  fc = fa;
  e = sb - sa;
  d = e;

  macheps = r8_epsilon ( );

  for ( ; ; )
  {
    if ( fabs ( fc ) < fabs ( fb ) )
    {
      sa = sb;
      sb = c;
      c = sa;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    tol = 2.0 * macheps * fabs ( sb ) + t;
    m = 0.5 * ( c - sb );

    if ( fabs ( m ) <= tol || fb == 0.0 )
    {
      break;
    }

    if ( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
    {
      e = m;
      d = e;
    }
    else
    {
      s = fb / fa;

      if ( sa == c )
      {
        p = 2.0 * m * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
        q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
      }

      if ( 0.0 < p )
      {
        q = - q;
      }
      else
      {
        p = - p;
      }

      s = e;
      e = d;

      if ( 2.0 * p < 3.0 * m * q - fabs ( tol * q ) &&
           p < fabs ( 0.5 * s * q ) )
      {
        d = p / q;
      }
      else
      {
        e = m;
        d = e;
      }
    }
    sa = sb;
    fa = fb;

    if ( tol < fabs ( d ) )
    {
      sb = sb + d;
    }
    else if ( 0.0 < m )
    {
      sb = sb + tol;
    }
    else
    {
      sb = sb - tol;
    }

    fb = cpp_temp_func(sb, y_ord, x_ord, n, n_grid, lambda2);

    if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
    {
      c = sa;
      fc = fa;
      e = sb - sa;
      d = e;
    }
  }
  return sb;
}


// [[Rcpp::export]]
double cpp_uniroot2(double lambda_min, double lambda_max,
                   arma::vec y_ord, arma::vec x_ord,
                   int n, int n_grid, double lambda2,
                   double tol = 1e-10, int max_iter = 300) {

  double b = lambda_max;
  double a = lambda_min;

  double fa = cpp_temp_func(a, y_ord, x_ord,
                            n, n_grid, lambda2);
  double fb = cpp_temp_func(b, y_ord, x_ord,
                            n, n_grid, lambda2);

  if(sign(fa) == 0) {
    return a;
  }
  if(sign(fb) == 0){
    return b;
  }

  if(sign(fa) == sign(fb)) {
    stop("The sign must be different at endpoints.");
  }

  int N = 1;
  double c = 0;

  while( N < max_iter) {
    c = (a + b)/2;
    double temp_func = cpp_temp_func(c, y_ord, x_ord,
                                     n, n_grid, lambda2);
    if(temp_func == 0 || (b-a)/2 < tol) {
      return c;
    }

    N = N + 1;
    if(sign(temp_func) == sign(fa)) {
      a = c;
      fa = temp_func;
    } else {
      b = c;
    }
  }

  return c;
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
                         int n, int n_grid,
                         double lam_tilde_old) {

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
  double lam_max = cpp_find_lamdaMax(y_ord, x_ord,
                                     n, n_grid, lam_tilde_old,
                                     lambda2);
  arma::vec tols(2);
  double tol = 1e-10;
  if(pow(lambda2,2) < tol) {
    tol = pow(lambda2,2);
  }

  if(lam_max < 0) {
    double b1 = as_scalar(cov(x_ord, y_ord))/as_scalar(var(x_ord));
    double b0 = as_scalar(mean(y_ord) - (b1*mean(x_ord)));
    f_hat = b0 + b1 * x_ord;
  } else {
    double temp_lam = cpp_uniroot3(0, lam_max, y_ord, x_ord,
                                  n, n_grid, lambda2,
                                  tol);
    f_hat = cpp_spline_raw(y_ord, x_ord, temp_lam, n, n_grid);
  }

  double temp_inner = as_scalar(1 - lambda1/(norm(f_hat)/pow(n, 0.5)) );
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

    arma::vec temp_ans = cpp_solve_prox(temp_y.elem(ord_mat.col(i)) - mean(temp_y),
                                        x_mat_ord.col(i), step_size*lambda1/n,
                                        step_size*lambda2/n, n, 500, 1);

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
List cpp_spp_one(arma::vec y, arma::mat x_ord,
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

  return List::create(Named("fhat") = new_ans(1),
                      Named("intercept") = as_scalar(new_ans(0)),
                      Named("conv") = converged );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

