#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"

/* Copyright: 2011 (dawei.zju@gmail.com) */
/* The Polytechnic University	version: 2011/09/07 */
/* Zhejiang University          version: 2012/02/27 */
/* The Polytechnic University	version: 2012/07/09 */
/* The Polytechnic University	version: 2012/07/11 */

/* Multiply and Accumulate a vector */
double vmac(double *a, double *b, int n) {
    double z = 0;
    int i=0;
    for(i=0;i<n;i++)
        z += a[i] * b[i];
    return z;
}

double sign(double num) {
    if( num < 0.0 )
        return -1.0F;
    else
        return +1.0F;
}

double vabs(double x) {
	if(x<0)
		return -1.0F * x;
	else
		return x;
}

double norm2(double *xr, double *xi, int n) {
    double xn = 0;
    int i;
    for (i=0;i<n;i++)
        xn += xr[i]*xr[i] + xi[i]*xi[i];
    return xn;
}

void ConjReverse(double *x, double *y, int N, double minus) {
    int i;
    for(i=0;i<N;i++)
        y[i] = x[N-1-i] * minus;
}

 
/* Rep: 1: M-CMA
		2: M-CMMA-1
		3: M-CMMA-16QAM
		4: M-CMMA-64QAM
		5: M-CMMA-4
		6: NCMA */
double errorfun(double a, double b, double *M, int N, int idx) {

	if (idx == 1)
		return M[0] - a*a;
	/* For 16-QAM */
    else if (idx == 2)
		return (1 + 9) - 2 * a*a * 9;
	/* For 16-QAM */
	/* For 64-QAM, reduce constellation to its proportional 16-QAM form */
    else if (idx == 3) {
		return 4.0F - vabs(a*a*9.0F-5.0F);
	}
	/* For 64-QAM, take 3 and 7 as two reference points, use same expression as #3 
	 * Reduce constellaion to its original points*/
	else if (idx == 4) {
        return 20 - vabs(a*a*49-29);
    }
	/* For 64-QAM, take 2 and 6 as two imaginary points, use same expression as #3 
	 * Reduce constellaion to its imaginary 16 points */
	else if (idx == 5) {
		return 16 - vabs(a*a*49-20);
    }
	else if (idx == 6) {
		double A1 = 4/7;
        double A2 = 2/7;
        double A3 = 1/7;
        double e1 = vabs(a) - A1;
        double e2 = vabs(e1) - A2;
        double e3 = vabs(e2) - A3;
        return -e3 * sign(e1) * sign(e2) * sign(a);
	}
	else if (idx == 7) {
        double e1 = 4/49 - (a*a-5/49);
		double e2 = 12/49 - (a*a-37/49);
		if (vabs(e1) < vabs(e2))
			return e1;
		else
			return e2;
	}
	else
		return 0.0F;
}

/* Updating of coefficients */
void updatecoeff(double *hr, double *hi, int Ntap,
        double *xr, double *xi, int Ndim,
        double yr, double yi, double mu, double *R, int NR, int IDX) {
    int i;
    double kr, ki;
    
    /* for these methods, the real and imag parts of y are multiplied by non-same err*/
    kr = mu * errorfun(yr, 0, R, NR, IDX) * yr;
    ki = mu * errorfun(yi, 0, R, NR, IDX) * yi;
    
    for(i=0;i<Ntap;i++) {
        *(hr+i) += kr * xr[i] + ki * xi[i];
        *(hi+i) += ki * xr[i] - kr * xi[i];
    }
}

/* Core CMA adaptive filter function */
void cmafilter(double *xr, double *xi, int Ndim,
        double *hr, double *hi, int Ntap,
        double mu, double *R, int NR,
        double *yr, double *yi, double *mse, double *det,
        double sps, bool dontskip, int id_error) {
    int i;
    int halftap = (Ntap-1)/2;
    int Ydim = Ndim-Ntap+1;
    double real=0, imag=0;
    
    for(i=0;i<Ydim;i++) {
       
        *(yr+i) = vmac(xr+i, hr, Ntap) - vmac(xi+i, hi, Ntap);
        *(yi+i) = vmac(xi+i, hr, Ntap) + vmac(xr+i, hi, Ntap);
        
        /* Updating filter coefficients */
        if (  dontskip || ( i % (int)sps == 0) ) {
            updatecoeff(hr,hi,Ntap,xr+i,xi+i,Ndim,yr[i],yi[i],mu,R,NR,id_error);
        }

        /* Calculate MSE and determinate of H */
        if (id_error == 2) {
			mse[i] = pow(yr[i],2) + pow(yi[i],2) - 2*R[0];
		}
        else {
			mse[i] = pow(errorfun(yr[i],yi[i],R,NR,id_error),1);
		}
                
        real = *(hr+halftap);
        imag = *(hi+halftap);      
        det[i] = sqrt(real*real + imag*imag);
        
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    double *xr			= mxGetPr(prhs[0]);             /* Input field vector, real part */
    double *xi			= mxGetPi(prhs[0]);             /* Input field vector, imag part */
    double *hr			= mxGetPr(prhs[1]);             /* 1st output filter, real part */
    double *hi			= mxGetPi(prhs[1]);             /* 1st output filter, imag part */
    int Ntap			= mxGetScalar(prhs[2]);         /* Length of filters, real scalar */
    double mu			= mxGetScalar(prhs[3]);         /* Alg. constant, real scalar */
    double *R			= mxGetPr(prhs[4]);             /* Convergence radii, real vector */
    double sps			= mxGetScalar(prhs[5]);         /* Samples x symbol, 1 or 2  */
    int id_error		= mxGetScalar(prhs[6]);			/* Error-function control flag */
    
    int Mdim			= mxGetM(prhs[0]);				/* Length of input field vector */
    int Npol			= mxGetN(prhs[0]);				/* Should be 2 */
    int Mfilter1		= mxGetM(prhs[1]);				/* Should be == Ntap */
    int Nfilter1		= mxGetN(prhs[1]);				/* Should be 2 */
    int NR				= mxGetM(prhs[4]);
    
    double *yr,*yi,*MSE,*DET;

    bool dontskip;
    
    if ( Ntap % 2 == 0)
        mexErrMsgTxt("Ntaps should be an ODD INTEGER.");
    
    if (!(Ntap == Mfilter1))
        mexErrMsgTxt("Filter length should be equal.");
    if (!(Npol == 1)||!(Nfilter1 == 1))
        mexErrMsgTxt("Single polarization is required.");
    
    /* Chechking the value of samples per symbol */
	dontskip = ( (int)sps == 1 );
    
    /* Chech if the vector that are supposed to be complex are really complex */
    /* If X, H1 or H2 are purely real numbers (which is often the case for H1 */
    /* and H2) then matlab does not allocates the memory for the imaginary    */
    /* part and the program will segfault if it tries to access it            */
    
    if (xi==NULL) {
        xi = mxCalloc(Mdim*Npol, sizeof(double));
        mxSetPi( (mxArray *)prhs[0], xi);
    }
    if (hi==NULL) {
        hi = mxCalloc(Ntap*Npol, sizeof(double));
        mxSetPi( (mxArray *)prhs[1], hi);
    }
    
    /* Next line allocates memory for a (Mdim-Ntap+1) x Npol complex vector */
    plhs[0] = mxCreateDoubleMatrix(Mdim-Ntap+1, Npol, mxCOMPLEX);
    yr = mxGetPr(plhs[0]);
    yi = mxGetPi(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(Mdim-Ntap+1, 1, mxREAL);
    MSE = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(Mdim-Ntap+1, 1, mxREAL);
    DET = mxGetPr(plhs[2]);
    
    cmafilter(xr,xi,Mdim,hr,hi,Ntap,mu,R,NR,yr,yi,MSE,DET,sps,dontskip,id_error);
    
    return;
}
