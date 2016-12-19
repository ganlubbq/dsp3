#include "mex.h"
#include "math.h"

/* Author: Wangdawei, 2010			*/
/* The Polytechnic University		*/

double makedecision(double x, double mn) {
    if ((int)mn == 4) {
        if (x >= 0)  return 1;
        else  return -1;
    }
    else if ((int)mn == 16) {
        if (x >= 0) {
            if (x >= 2) return 3;
            else return 1;  }
        else {
            if (x <= -2) return -3;
            else return -1;  }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,  const mxArray *prhs[]) {
    double *xr  = mxGetPr(prhs[0]);     /* Input field vector, real part */
    double *xi  = mxGetPi(prhs[0]);     /* Input field vector, imag part */
    double mu = mxGetScalar(prhs[1]);   /* Step parameter  */
    double *phi = mxGetPr(prhs[2]);
    double mn = mxGetScalar(prhs[3]);
    
    int M     = mxGetM(prhs[0]);     /* Length of input field vector  */
    int N     = mxGetN(prhs[0]);     /* Should be 2                   */
    
    double *yr;
    double *yi;
    
    int i,k;
    double rerr, ierr;
    
    if (xi==NULL) {   /* If necessary allocates memory for the imaginary part of X */
        xi = mxCalloc(M*N, sizeof(double));
        mxSetPi( (mxArray *)prhs[0], xi);
    }
   
    plhs[0] = mxCreateDoubleMatrix(M, N, mxCOMPLEX);
    yr = mxGetPr(plhs[0]);     /* Output field vector, real part */
    yi = mxGetPi(plhs[0]);     /* Output field vector, imag part */
    for (k = 0;k < N; k++) {
        for (i = 0;i < M; i++) {
            yr[i+k*M] = xr[i+k*M]*cos(phi[k]) + xi[i+k*M]*sin(phi[k]);
            yi[i+k*M] = xi[i+k*M]*cos(phi[k]) - xr[i+k*M]*sin(phi[k]);
            rerr = yr[i+k*M] - makedecision(yr[i+k*M], mn);
            ierr = yi[i+k*M] - makedecision(yi[i+k*M], mn);
            phi[k] -= mu * (yi[i+k*M]*rerr - yr[i+k*M]*ierr);
        }
    }
    return;
}
