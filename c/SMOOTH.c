#include "mex.h"

/* Copyright: 2011 (dawei.zju@gmail.com) */
/* The Polytechnic University	version: 2012/08/29 */

void vrot(double *yr, double *yi, double *xr, double *xi, int n, double p) {
    int i = 0;
    for (i=0;i<n;i++) {
        yr[i] += xr[i]*cos(p) - xi[i]*sin(p);
        yi[i] += xr[i]*sin(p) + xi[i]*cos(p);
    }
}

/* bs must be odd */
void smooth_complex(double *yr, double *yi, double *xr, double *xi, int Mlen, int Npol, int bs) {
    int half_bs,k,kk,N,i;
    N = Mlen*Npol;
    half_bs = (bs-1)/2;

    for (k=0;k<N;k++) {
        yr[k] = 0;
        yi[k] = 0;
        if (k<half_bs) {
            for (kk=0;kk<2*k+1;kk++) {
                yr[k] += xr[kk];
                yi[k] += xi[kk];
            }
            yr[k] = yr[k] / (2*k+1);
            yi[k] = yi[k] / (2*k+1);
        }
        else if (k>(N-half_bs-1)) {
            for (kk=2*k+1-N;kk<N;kk++) {
                yr[k] += xr[kk];
                yi[k] += xi[kk];
            }
            yr[k] = yr[k] / (2*N-2*k-1);
            yi[k] = yi[k] / (2*N-2*k-1);
        }
        else {
            for (kk=k-half_bs;kk<k+half_bs+1;kk++) {
                yr[k] += xr[kk];
                yi[k] += xi[kk];
            }
            yr[k] = yr[k] / bs;
            yi[k] = yi[k] / bs;
        }
        
    }

}
/* bs must be odd */
void smooth_real(double *yr, double *xr, int Mlen, int Npol, int bs) {
    int half_bs,k,kk,N;
    N = Mlen*Npol;
    half_bs = (bs-1)/2;
    
    for (k=0;k<N;k++) {
        yr[k] = 0;
        if (k<half_bs) {
            for (kk=0;kk<2*k+1;kk++) {
                yr[k] += xr[kk];
            }
            yr[k] = yr[k] / (2*k+1);
        }
        else if (k>(N-half_bs-1)) {
            for (kk=2*k+1-N;kk<N;kk++) {
                yr[k] += xr[kk];
            }
            yr[k] = yr[k] / (2*N-2*k-1);
        }
        else {
            for (kk=k-half_bs;kk<k+half_bs+1;kk++) {
                yr[k] += xr[kk];
            }
            yr[k] = yr[k] / bs;
        }
        
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *xr = mxGetPr(prhs[0]);          /* Input field vector, real part */
    double *xi = mxGetPi(prhs[0]);          /* Input field vector, imag part */
    int bs = mxGetScalar(prhs[1]);
    
    
    int Mlen = mxGetM(prhs[0]);             /* Length of input field vector  */
    int Npol = mxGetN(prhs[0]);             /* Should be 2                   */
    
    int i;
    double *yr,*yi;

    
    if (bs%2 == 0)
        bs += 1;
   
  
    if (xi==NULL) {
        plhs[0] = mxCreateDoubleMatrix(Mlen, Npol, mxREAL);
        yr = mxGetPr(plhs[0]);
        smooth_real(yr, xr, Mlen, Npol, bs);
    }
    else {
        plhs[0] = mxCreateDoubleMatrix(Mlen, Npol, mxCOMPLEX);
        yr = mxGetPr(plhs[0]);
        yi = mxGetPi(plhs[0]);
        smooth_complex(yr, yi, xr, xi, Mlen, Npol, bs);
    }
    return;
}
