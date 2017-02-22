#include<stdlib.h>
#include<iostream>
#include<math.h>

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

void callSS_Fortran(double* y, double* x, double *sy, double *lambda, int *n_point)
{
  double lam = (*lambda);
  double n = (*n_point);
  double h = log(2*lam/n); //log lambda
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
  int ngrid = 0;
  double xg = 0;
  double yg = 0;

  int job[3]; //Assumes Xs are sorted
  job[0] = 0; job[1] = 0; job[2] = 1;
  int ideriv = 0;
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
       &ngrid,
       &xg,
       &yg,
       job,
       &ideriv,
       &ierr);
}
}
