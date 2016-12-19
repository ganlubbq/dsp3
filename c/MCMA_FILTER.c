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
    int i;
    for(i=0;i<n;i++)
        z += a[i] * b[i];
    return z;
}

double vabs(double x) {
	if(x<0)
		return -1.0F * x;
	else
		return x;
}

double sign(double num) {
    if( num < 0.0 )
        return -1.0F;
    else
        return +1.0F;
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
 	*	2: M-CMMA-1
	*	3: M-CMMA-16QAM
	*	4: M-CMMA-64QAM
	*	5: M-CMMA-4
	*	6: NCMA */
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
        double *xr, double *xi, int Mdim,
        double yr, double yi, double mu, double *R, int NR, int IDX) {
    int i;
    double kr, ki, x_norm;
 
    kr = mu * errorfun(yr, 0, R, NR, IDX) * yr;
    ki = mu * errorfun(yi, 0, R, NR, IDX) * yi;
    
    for( i=0;i<Ntap;i++) {
        *(hr+i) += kr * xr[i] + ki * xi[i];
        *(hi+i) += ki * xr[i] - kr * xi[i];
        *(hr+i+Ntap) += kr * xr[i+Mdim] + ki * xi[i+Mdim];
        *(hi+i+Ntap) += ki * xr[i+Mdim] - kr * xi[i+Mdim];
    }
}

/* Core CMA adaptive filter function */
void cmafilter(double *xr, double *xi, int Mdim,
        double *h1r, double *h1i, double *h2r, double *h2i, int Ntap,
        double mu, double *R, int NR,
        double *yr, double *yi, double *mse, double *det,
        double sps, bool dontskip,
        int id_error, int id_stage, int id_method) {

    int i, j, halftap, Ydim;
    double real, imag, msex, msey;
    double *h1rv, *h1iv, *h2rv, *h2iv;
    
    halftap = (Ntap-1)/2;
    Ydim = Mdim-Ntap+1;
    h1rv = mxMalloc(2*Ntap*sizeof(double));
    h2rv = mxMalloc(2*Ntap*sizeof(double));
    h1iv = mxMalloc(2*Ntap*sizeof(double));
    h2iv = mxMalloc(2*Ntap*sizeof(double));
    
    for( i=0;i<Ydim;i++) {
       
        /* Output first column real part yr = Real( x * h1) */
        *(yr+i) = vmac(xr+i, h1r, Ntap) - vmac(xi+i, h1i, Ntap)             /* 1st filter column */
        +vmac(xr+i+Mdim, h1r+Ntap, Ntap)-vmac(xi+i+Mdim, h1i+Ntap, Ntap);   /* 2nd filter column */;
        /* Output first column imag part yi = Imag( x * h1) */
        *(yi+i) = vmac(xi+i, h1r, Ntap) + vmac(xr+i, h1i, Ntap)             /* 1st filter column */
        +vmac(xi+i+Mdim, h1r+Ntap, Ntap)+vmac(xr+i+Mdim, h1i+Ntap, Ntap);   /* 2nd filter column */;
        /* Output second column real part yr = Real( x * h2) */
        *(yr+i+Ydim) = vmac(xr+i, h2r, Ntap) - vmac(xi+i, h2i, Ntap)        /* 1st filter column */
        +vmac(xr+i+Mdim, h2r+Ntap, Ntap)-vmac(xi+i+Mdim, h2i+Ntap, Ntap);   /* 2nd filter column */;
        /* Output second column imag part yi = Imag( x * h1) */
        *(yi+i+Ydim) = vmac(xi+i, h2r, Ntap) + vmac(xr+i, h2i, Ntap)        /* 1st filter column */
        +vmac(xi+i+Mdim, h2r+Ntap, Ntap)+vmac(xr+i+Mdim, h2i+Ntap, Ntap);   /* 2nd filter column */;
        
        /* Updating filter coefficients */
        if ( id_stage > 1) {
            if (  dontskip || ( i % (int)sps == 0)) {
                updatecoeff(h1r,h1i,Ntap,xr+i,xi+i,Mdim,yr[i],yi[i],mu,R,NR,id_error);
                updatecoeff(h2r,h2i,Ntap,xr+i,xi+i,Mdim,yr[i+Ydim],yi[i+Ydim],mu,R,NR,id_error);
            }
        }
        else {
            if (  dontskip || ( i % (int)sps == 0)) {
                updatecoeff(h1r, h1i,Ntap,xr+i,xi+i,Mdim,yr[i],yi[i],mu,R,NR,id_error);
                ConjReverse(h1r, h1rv, Ntap, +1.0F); 
                ConjReverse(h1r+Ntap, h1rv+Ntap, Ntap, +1.0F);
                ConjReverse(h1i, h1iv, Ntap, -1.0F); 
                ConjReverse(h1i+Ntap, h1iv+Ntap, Ntap, -1.0F);
                for(j=0;j<Ntap;j++) {
                    *(h2r+Ntap+j) = *(h1rv+j);
                    *(h2i+Ntap+j) = *(h1iv+j);
                    *(h2r+j) = *(h1rv+Ntap+j) * -1.0F;
                    *(h2i+j) = *(h1iv+Ntap+j) * -1.0F;
                }
            }
        }

        /* Calculate error signal, upon convergence, error signal should approach zero */
		msex = errorfun(yr[i],0,R,NR,id_error)*yr[i] + errorfun(yi[i],0,R,NR,id_error)*yi[i];
		msey = errorfun(yr[i+Ydim],0,R,NR,id_error)*yr[i+Ydim] + errorfun(yi[i+Ydim],0,R,NR,id_error)*yi[i+Ydim];
		mse[i] = (msex+msey)/2;

		/* Calculate MSE and determinate of H */
		/*msex = pow(yr[i],2)-R[0] + pow(yi[i],2)-R[0];
		msey = pow(yr[i+Ydim],2)-R[0] + pow(yi[i+Ydim],2)-R[0];
        mse[i] = pow(msex,2) + pow(msey,2);*/
                
        real = (*(h1r+halftap)* *(h2r+Ntap+halftap)-*(h1i+halftap)* *(h2i+Ntap+halftap))
			-(*(h1r+Ntap+halftap)* *(h2r+halftap)-*(h1i+Ntap+halftap)* *(h2i+halftap));
        imag = (*(h1r+halftap)* *(h2i+Ntap+halftap)+*(h1i+halftap)* *(h2r+Ntap+halftap))
			-(*(h1r+Ntap+halftap)* *(h2i+halftap)+*(h1i+Ntap+halftap)* *(h2r+halftap));        
        det[i] = sqrt(real*real + imag*imag);
        
        /* Anti-sigularity */
        if ( id_method == 1) {
            if ( *(det+i) < 0.01) {
                ConjReverse(h1r, h1rv, Ntap, +1.0F); 
                ConjReverse(h1r+Ntap, h1rv+Ntap, Ntap, +1.0F);
                ConjReverse(h1i, h1iv, Ntap, -1.0F); 
                ConjReverse(h1i+Ntap, h1iv+Ntap, Ntap, -1.0F);
                for(j=0;j<Ntap;j++) {
                    *(h2r+Ntap+j) = *(h1rv+j);
                    *(h2i+Ntap+j) = *(h1iv+j);
                    *(h2r+j) = *(h1rv+Ntap+j) * -1.0F;
                    *(h2i+j) = *(h1iv+Ntap+j) * -1.0F;
                }
            }
        }
        if ( id_method == 2) {
            if ( !((i+1) % 100)) {
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
    
    double *xr		= mxGetPr(prhs[0]);             /* Input field vector, real part */
    double *xi		= mxGetPi(prhs[0]);             /* Input field vector, imag part */
    double *h1r		= mxGetPr(prhs[1]);             /* 1st output filter, real part */
    double *h1i		= mxGetPi(prhs[1]);             /* 1st output filter, imag part */
    double *h2r		= mxGetPr(prhs[2]);             /* 2nd output filter, real part */
    double *h2i		= mxGetPi(prhs[2]);             /* 2nd output filter, imag part */
    int Ntap		= mxGetScalar(prhs[3]);         /* Length of filters, real scalar */
    double mu		= mxGetScalar(prhs[4]);         /* Alg. constant, real scalar */
    double *R		= mxGetPr(prhs[5]);             /* Convergence radii, real vector */
    double sps		= mxGetScalar(prhs[6]);         /* Samples x symbol, 1 or 2 */
    int id_error	= mxGetScalar(prhs[7]);			/* Error function control flag */
    int id_stage	= mxGetScalar(prhs[8]);			/* Stage control flag */
    int id_method	= mxGetScalar(prhs[9]);			/* Method control flag */
    
    int Mdim		= mxGetM(prhs[0]);				/* Length of input field vector */
    int Npol		= mxGetN(prhs[0]);				/* Should be 2 */
    int Mfilt_1		= mxGetM(prhs[1]);				/* Should be == Ntap */
    int Nfilt_1		= mxGetN(prhs[1]);				/* Should be 2 */
    int Mfilt_2		= mxGetM(prhs[2]);				/* Should be == Ntap */
    int Nfilt_2		= mxGetN(prhs[2]);				/* Should be 2 */
    int NR			= mxGetM(prhs[5]);

	double *yr,*yi,*MSE,*DET;

    bool dontskip;

    if ( Ntap % 2 == 0)
        mexErrMsgTxt("Ntaps should be an ODD INTEGER.");
    
    if ( id_error<1 || id_error>7)
        mexErrMsgTxt("Error function ID should be 1~7");
    
    if ( !(Ntap == Mfilt_1) || !(Ntap == Mfilt_2))
        mexErrMsgTxt("Filter length should be equal.");
    
    if ( !(Npol == 2) || !(Nfilt_1 == 2) || !(Nfilt_2 == 2))
        mexErrMsgTxt("Dual polarizations are required.");
    
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
    
    cmafilter(xr,xi,Mdim,
            h1r,h1i,h2r,h2i,Ntap,
            mu,R,NR,
            yr,yi,MSE,DET,
            sps,dontskip,
            id_error,id_stage,id_method);
    
    return;
}
