#include "mex.h"
#include <cmath>

using namespace std;

extern "C"{
  //blas functions
  double dgemv_(char* TRANS, const int* M, const int* N, 
		double alpha, double* A, const int* LDA, 
		double* X, const int* INCX, double* beta,
 		double* Y, const int* INCY);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //#define A_in		prhs[0]
  #define X_in		prhs[0]
  #define Y_out		plhs[0]

  int m = mxGetM(A_in);
  int n = mxGetN(X_in);
  char trans = 'N';
  if (nlhs == 2) trans = 'T';
  const int incx = 1, incy = 1;
  double alpha = 1.
  double beta= 0.;

  Y_out = mxCreateDoubleMatrix(m, n, mxREAL);
  Y = mxGetPr(y_out);
  // A = 
  X = mxGetPr(x_in);
  
  dgemv_(&trans, &m, &n, &alpha, A, &n, X, 
	&incx, &beta, Y, &incy);
  return;

}
