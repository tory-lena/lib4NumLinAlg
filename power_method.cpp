#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

using namespace std;

// FORTRAN adds _ after all the function names
// and all variables are called by reference
extern "C"{
  // BLAS functions
  double dscal_(const int* N, double* alpha, double* X, const int* INCX);

  double dcopy_(const int* N, double* X, const int* INCX, double* Y, const int* INCY);	
  double ddot_(const int* N, double* X, const int* INCX, double* Y, const int* INCY);

  double dgemv_(char* TRANS, const int* M, const int* N, double* alpha,
               double* A, const int* LDA, double* X, const int* INCX,
               double* beta, double* Y, const int* INCY);
  
  double dnrm2_(const int* N, const double *X, const int* INCX);
}

void read_matrix(ifstream& file, int n, double * A) {
  for (int i = 0; i < n; i++ ) {
    for (int j = 0; j < n; j++) {
      // Write in column major order
      file >> A[j*n+i];
    }
  }
}

int main( int argc, char** argv ){
  if (argc != 2) {
    printf("usage: power_method matrix\n");
    printf("  matrix should be text file of form:\n");
    printf("    n\n");
    printf("    n lines of n space-delimited numbers\n");

    return 0;
  }

  // Get matrix
  ifstream mat_file( argv[1] );
  if (!mat_file) {
    printf("Could not open file.");
    return 1;
  }

  int n;
  mat_file >> n;
  double *A = (double *) malloc(n * n * sizeof(double));
  read_matrix(mat_file, n, A);

  // Random number generator
  random_device rd;
  mt19937 gen(rd());
  std::uniform_real_distribution<double> unif(-1., 1.);

  // implement power method here
  double lambda_old = 0;
  double lambda = 1;
  int k = 1;
  char trans = 'N';
  double tol = 1E-6;
  int max_it = 20;
  const int incx = 1, incy = 1;
  //double dot, norm;
  double alpha, beta= 0.;

  double *x = (double *) calloc(n, sizeof(double));
  for (int i = 0; i<n; i++){
    x[i] = unif(gen);
  }
  double *y = (double *) calloc(n, sizeof(double));

  while ((abs(lambda - lambda_old) > tol) && (k < max_it)){

    lambda_old = lambda;
    alpha = 1./ddot_(&n, x, &incx, x, &incy);

    dgemv_(&trans, &n, &n, &alpha, A, &n, x, &incx, &beta, y, &incy);
    lambda = ddot_(&n, x, &incx, y, &incy);
    // y <- alpha*A*x + beta*y
    alpha = 1.;
    dgemv_(&trans, &n, &n, &alpha, A, &n, x, &incx, &beta, y, &incy);
    alpha = 1./dnrm2_(&n, y, &incx);
    dscal_(&n, &alpha, y, &incx);
    dcopy_(&n, y, &incx, x, &incy);
    k++;
  }
  cout << abs(lambda - lambda_old) << endl;
  cout << k << endl;

  // Print results
  printf("Eigenvalue : %f\n", lambda);
  printf("Eigenvector:\n");
  for (int i = 0; i < n; i++)
    printf("%f\n", x[i]);;

  // cleanup
  free(x);

  return 0;


};
