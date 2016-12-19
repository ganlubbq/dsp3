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

 
/* Rep: 1.CMA. 2.MCMA. 3.CMMA. 4.RD-CMA. 5&6.New CMA. 7 NCMA */
double errorfun(double a, double b, double *M, int N, double idx) {

    if ((int)idx == 1 || (int)idx == 7) 
        return M[0] - a*a - b*b;
    
	else if ((int)idx == 2)
		return M[0] - a*a;
    
    else if ((int)idx == 3) {
        double A1 = 0.5*(M[0]+M[1]);
        double A2 = 0.5*(M[2]-M[0]);
        double A3 = 0.5*(M[2]-M[1]);
        double e1 = sqrt(a*a+b*b)-A1;
        double e2 = vabs(e1)-A2;
        double e3 = vabs(e2)-A3;
        return -e3 * sign(e1) * sign(e2) / sqrt(a*a+b*b);
    }
    else if ((int)idx == 4) {
        int k;
        double mini=10000, tmp;
        for (k=0;k<N;k++) {
            tmp = M[k]-a*a-b*b;
            if ( vabs(tmp) < vabs(mini) )
                    mini = tmp;
        }
        return mini;
    }
    else if ((int)idx == 5)
		return (1*1 + 3*3) - 2 * a*a * 9;
    
    else if ((int)idx == 6)
		return 4 - vabs(a*a*9-5);
    
	else if ((int)idx == 8) {
		double A1,A2,A3,e1,e2,e[3];
		double mini = 10000;
        int k;
        A1 = 0.5*(M[0]+M[1]);
        A2 = 0.5*(M[2]-M[0]);
        A3 = 0.5*(M[2]-M[1]);
		e1 = sqrt(a*a+b*b)-A1;
		e2 = vabs(e1)-A2;
        e[0] = -(vabs(e2)-A3) * sign(e1) * sign(e2) / sqrt(a*a+b*b);
		A1 = 0.5*(M[3]+M[4]);
        A2 = 0.5*(M[5]-M[3]);
        A3 = 0.5*(M[5]-M[4]);
        e1 = sqrt(a*a+b*b)-A1;
		e2 = vabs(e1)-A2;
        e[1] = -(vabs(e2)-A3) * sign(e1) * sign(e2) / sqrt(a*a+b*b);
		A1 = 0.5*(M[6]+M[7]);
        A2 = 0.5*(M[8]-M[6]);
        A3 = 0.5*(M[8]-M[7]);
        e1 = sqrt(a*a+b*b)-A1;
		e2 = vabs(e1)-A2;
        e[2] = -(vabs(e2)-A3) * sign(e1) * sign(e2) / sqrt(a*a+b*b);
		for (k=0;k<3;k++) {
            if ( vabs(e[k]) < vabs(mini) )
                    mini = e[k];
        }
        return mini;
    }
	else if ((int)idx == 9) {
        double A1 = 4/7;
        double A2 = 2/7;
        double A3 = 1/7;
        double e1 = vabs(a) - A1;
        double e2 = vabs(e1) - A2;
        double e3 = vabs(e2) - A3;
        return -e3 * sign(e1) * sign(e2) * sign(a);
    }
	else
		return 0.0F;
}

/* Updating of coefficients */
void updatecoeff(double *hr, double *hi, int Ntap,
        double *xr, double *xi, int Ndim,
        double yr, double yi, double mu, double *R, int NR, double IDX) {
    int i;
    double kr, ki, x_norm;
    
    /* for these methods, the real and imag parts of y are multiplied by same err */
    if ((int)IDX == 1 || (int)IDX == 3 || (int)IDX == 4 || (int)IDX == 8) {
        kr = mu * errorfun(yr, yi, R, NR, IDX) * yr;
        ki = mu * errorfun(yr, yi, R, NR, IDX) * yi;
    }
    
    /* for these methods, the real and imag parts of y are multiplied by non-same err*/
    else if ((int)IDX == 2 || (int)IDX == 5 || (int)IDX == 6 || (int)IDX == 9) {
        kr = mu * errorfun(yr, 0, R, NR, IDX) * yr;
        ki = mu * errorfun(yi, 0, R, NR, IDX) * yi;
    }

    else if ((int)IDX == 7) {
        x_norm = norm2(xr,xi,Ntap);
        kr = mu * errorfun(yr, yi, R, NR, IDX) * yr / x_norm;
        ki = mu * errorfun(yr, yi, R, NR, IDX) * yi / x_norm;
    }
    
    for(i=0;i<Ntap;i++) {
        *(hr+i) += kr * xr[i] + ki * xi[i];
        *(hi+i) += ki * xr[i] - kr * xi[i];
        *(hr+i+Ntap) += kr * xr[i+Ndim] + ki * xi[i+Ndim];
        *(hi+i+Ntap) += ki * xr[i+Ndim] - kr * xi[i+Ndim];
    }
}

/* Core CMA adaptive filter function */
void cmafilter(double *xr, double *xi, int Ndim,
        double *h1r, double *h1i, double *h2r, double *h2i, int Ntap,
        double mu, double *R, int NR,
        double *yr, double *yi, double *mse, double *det,
        double sps, bool dontskip,
        double id_error, double id_stage, double id_method) {
    int i, j;
    int halftap = (Ntap-1)/2;
    int Ydim = Ndim-Ntap+1;
    double real=0, imag=0, msex=0, msey=0;
    double *h1rv, *h1iv, *h2rv, *h2iv;

    h1rv = mxMalloc(2*Ntap*sizeof(double));
    h2rv = mxMalloc(2*Ntap*sizeof(double));
    h1iv = mxMalloc(2*Ntap*sizeof(double));
    h2iv = mxMalloc(2*Ntap*sizeof(double));
    
    for(i=0;i<Ydim;i++) {
       
        /* Output first column real part yr = Real( x * h1) */
        *(yr+i) = vmac(xr+i, h1r, Ntap) - vmac(xi+i, h1i, Ntap)             /* 1st filter column */
        +vmac(xr+i+Ndim, h1r+Ntap, Ntap)-vmac(xi+i+Ndim, h1i+Ntap, Ntap);   /* 2nd filter column */;
        /* Output first column imag part yi = Imag( x * h1) */
        *(yi+i) = vmac(xi+i, h1r, Ntap) + vmac(xr+i, h1i, Ntap)             /* 1st filter column */
        +vmac(xi+i+Ndim, h1r+Ntap, Ntap)+vmac(xr+i+Ndim, h1i+Ntap, Ntap);   /* 2nd filter column */;
        /* Output second column real part yr = Real( x * h2) */
        *(yr+i+Ydim) = vmac(xr+i, h2r, Ntap) - vmac(xi+i, h2i, Ntap)        /* 1st filter column */
        +vmac(xr+i+Ndim, h2r+Ntap, Ntap)-vmac(xi+i+Ndim, h2i+Ntap, Ntap);   /* 2nd filter column */;
        /* Output second column imag part yi = Imag( x * h1) */
        *(yi+i+Ydim) = vmac(xi+i, h2r, Ntap) + vmac(xr+i, h2i, Ntap)        /* 1st filter column */
        +vmac(xi+i+Ndim, h2r+Ntap, Ntap)+vmac(xr+i+Ndim, h2i+Ntap, Ntap);   /* 2nd filter column */;
        
        /* Updating filter coefficients */
        if ((int)id_stage > 1) {
            if (  dontskip || ( i % (int)sps == 0) ) {
                updatecoeff(h1r, h1i,Ntap,xr+i,xi+i,Ndim,yr[i],     yi[i],     mu,R,NR,id_error);
                updatecoeff(h2r, h2i,Ntap,xr+i,xi+i,Ndim,yr[i+Ydim],yi[i+Ydim],mu,R,NR,id_error);
            }
        }
        else {
            if (  dontskip || ( i % (int)sps == 0) ) {
                
                updatecoeff(h1r, h1i,Ntap,xr+i,xi+i,Ndim,yr[i],yi[i],mu,R,NR,id_error);
                
                ConjReverse(h1r,      h1rv,      Ntap, +1.0F); 
                ConjReverse(h1r+Ntap, h1rv+Ntap, Ntap, +1.0F);
                ConjReverse(h1i,      h1iv,      Ntap, -1.0F); 
                ConjReverse(h1i+Ntap, h1iv+Ntap, Ntap, -1.0F);
                
                for(j=0;j<Ntap;j++) {
                    *(h2r+j+Ntap) = *(h1rv+j);
                    *(h2i+j+Ntap) = *(h1iv+j);
                    *(h2r+j)      = *(h1rv+j+Ntap) * -1.0F;
                    *(h2i+j)      = *(h1iv+j+Ntap) * -1.0F;
                }
            }
        }

        /* Calculate MSE and determinate of H */
        if ((int)id_error == 2) {
			msex = pow(yr[i],     2)-R[0] + pow(yi[i],     2)-R[0];
			msey = pow(yr[i+Ydim],2)-R[0] + pow(yi[i+Ydim],2)-R[0];
		}
        else {
			msex = errorfun(yr[i],     yi[i],     R, NR, id_error);
			msey = errorfun(yr[i+Ydim],yi[i+Ydim],R, NR, id_error);
		}
        mse[i] = pow(msex,2) + pow(msey,2);
                
        real = (*(h1r+halftap)* *(h2r+Ntap+halftap)-*(h1i+halftap)* *(h2i+Ntap+halftap))
			-(*(h1r+Ntap+halftap)* *(h2r+halftap)-*(h1i+Ntap+halftap)* *(h2i+halftap));
        imag = (*(h1r+halftap)* *(h2i+Ntap+halftap)+*(h1i+halftap)* *(h2r+Ntap+halftap))
			-(*(h1r+Ntap+halftap)* *(h2i+halftap)+*(h1i+Ntap+halftap)* *(h2r+halftap));
        
        det[i] = sqrt(real*real + imag*imag);
        
        /* anti-sigularity */
        if ( id_method == 1 ) {
            if ( *(det+i) < 0.01 ) {
                ConjReverse(h1r,      h1rv,      Ntap, +1.0F); 
                ConjReverse(h1r+Ntap, h1rv+Ntap, Ntap, +1.0F);
                ConjReverse(h1i,      h1iv,      Ntap, -1.0F); 
                ConjReverse(h1i+Ntap, h1iv+Ntap, Ntap, -1.0F);
                for(j=0;j<Ntap;j++) {
                    *(h2r+j+Ntap) = *(h1rv+j);
                    *(h2i+j+Ntap) = *(h1iv+j);
                    *(h2r+j)      = *(h1rv+j+Ntap) * -1.0F;
                    *(h2i+j)      = *(h1iv+j+Ntap) * -1.0F;
                }
            }
        }
        if ( id_method == 2 ) {
            if ( !((i+1) % 100) ) {
                ConjReverse(h1r, h1rv, Ntap, +1.0F); 
                ConjReverse(h1r+Ntap, h1rv+Ntap, Ntap, +1.0F);
                ConjReverse(h2r, h2rv, Ntap, +1.0F); 
                ConjReverse(h2r+Ntap, h2rv+Ntap, Ntap, +1.0F);
                ConjReverse(h1i, h1iv, Ntap, -1.0F); 
                ConjReverse(h1i+Ntap, h1iv+Ntap, Ntap, -1.0F);
                ConjReverse(h2i, h2iv, Ntap, -1.0F); 
                ConjReverse(h2i+Ntap, h2iv+Ntap, Ntap, -1.0F);
                for(j=0;j<Ntap;j++) {
                    *(h1r+j) = 0.5*(*(h1r+j)+*(h2rv+Ntap+j));
                    *(h1i+j) = 0.5*(*(h1i+j)+*(h2iv+Ntap+j));
                    *(h2r+Ntap+j) = 0.5*(*(h2r+Ntap+j)+*(h1rv+j));
                    *(h2i+Ntap+j) = 0.5*(*(h2i+Ntap+j)+*(h1iv+j));
                    *(h1r+Ntap+j) = 0.5*(*(h1r+Ntap+j)-*(h2rv+j));
                    *(h1i+Ntap+j) = 0.5*(*(h1i+Ntap+j)-*(h2iv+j));
                    *(h2r+j) = 0.5*(*(h2r+j)-*(h1rv+Ntap+j));
                    *(h2i+j) = 0.5*(*(h2i+j)-*(h1iv+Ntap+j));
                }
            }
        }
        
    }
    mxFree(h1rv);
	mxFree(h2rv);
	mxFree(h1iv);
	mxFree(h2iv);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    double *xr			= mxGetPr(prhs[0]);             /* Input field vector, real part */
    double *xi			= mxGetPi(prhs[0]);             /* Input field vector, imag part */
    double *h1r			= mxGetPr(prhs[1]);             /* 1st output filter, real part */
    double *h1i			= mxGetPi(prhs[1]);             /* 1st output filter, imag part */
    double *h2r			= mxGetPr(prhs[2]);             /* 2nd output filter, real part */
    double *h2i			= mxGetPi(prhs[2]);             /* 2nd output filter, imag part */
    int Ntap			= mxGetScalar(prhs[3]);         /* Length of filters, real scalar */
    double mu			= mxGetScalar(prhs[4]);         /* Alg. constant, real scalar */
    double *R			= mxGetPr(prhs[5]);             /* Convergence radii, real vector */
    double sps			= mxGetScalar(prhs[6]);         /* Samples x symbol, 1 or 2  */
    double id_error		= mxGetScalar(prhs[7]);			/* Error function control flag */
    double id_stage		= mxGetScalar(prhs[8]);			/* Stage control flag  */
    double id_method	= mxGetScalar(prhs[9]);			/* Method control flag  */
    
    int Mdim			= mxGetM(prhs[0]);				/* Length of input field vector */
    int Npol			= mxGetN(prhs[0]);				/* Should be 2 */
    int Mfilter1		= mxGetM(prhs[1]);				/* Should be == Ntap */
    int Nfilter1		= mxGetN(prhs[1]);				/* Should be 2 */
    int Mfilter2		= mxGetM(prhs[2]);				/* Should be == Ntap */
    int Nfilter2		= mxGetN(prhs[2]);				/* Should be 2 */
    int NR				= mxGetM(prhs[5]);
    
    double *yr,*yi,*MSE,*DET;
    
    bool dontskip;
    
    if (Ntap % 2 == 0)
        mexErrMsgTxt("Ntaps should be an ODD INTEGER.");
    
    if (id_error<1 || id_error>9)
        mexErrMsgTxt("error function ID should be 1~9");
    
    if (!(Ntap == Mfilter1)||!(Ntap == Mfilter2))
        mexErrMsgTxt("Filter length should be equal.");
    
    if (!(Npol == 2)||!(Nfilter1 == 2)||!(Nfilter2 == 2))
        mexErrMsgTxt("Two polarizations needed.");
    
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
    if (h1i==NULL) {
        h1i = mxCalloc(Ntap*Npol, sizeof(double));
        mxSetPi( (mxArray *)prhs[1], h1i);
    }
    if (h2i==NULL) {
        h2i = mxCalloc(Ntap*Npol, sizeof(double));
        mxSetPi( (mxArray *)prhs[2], h2i);
    }
    
    /* Next line allocates memory for a (Mdim-Ntap+1) x Npol complex vector */
    plhs[0] = mxCreateDoubleMatrix(Mdim-Ntap+1, Npol, mxCOMPLEX);
    yr = mxGetPr(plhs[0]);
    yi = mxGetPi(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(Mdim-Ntap+1, 1, mxREAL);
    MSE = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(Mdim-Ntap+1, 1, mxREAL);
    DET = mxGetPr(plhs[2]);
    
    cmafilter(xr,xi,Mdim,h1r,h1i,h2r,h2i,Ntap,mu,R,NR,yr,yi,MSE,DET,sps,dontskip,id_error,id_stage,id_method);
    
    return;
}
