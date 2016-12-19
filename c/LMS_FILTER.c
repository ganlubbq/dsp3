#include <math.h>
#include "mex.h"

/* Author: Wangdawei, 2010			*/
/* The Polytechnic University	version: 2011/09/07 office pc	*/

/* Multiply and Accumulate a vector */
double vmac(double *a, double *b, int n)
{
    double z = 0;
    int i=0;
    for(i=0;i<n;i++)
        z += a[i] * b[i];
    return z;
}
/* Fast square root, from Quake III */
double Q_sqrt(float num)
{
    long i;
    float x2,y;
    const float threehalfs = 1.5F;
    
    x2 = num*0.5F;
    y = num;
    i = *(long*)&y;
    i = 0x5f3759df - (i>>1);
    y = *(float*)&i;
    y = y*(threehalfs-(x2*y*y));
    y = y*(threehalfs-(x2*y*y));
    
    return 1.0F/y;
}
double Q_abs(double num)
{
    if(num<0)
        return -1.0F * num;
    else
        return num;
}



void errorfun(double a,double b, double *Mr, double *Mi, int N, double *err)
{
    double min = 10000000, distance = 0;
    int kkk = 0, i;
    for ( i=0; i<N; i++) {
        distance = Q_sqrt((Mr[i]-a)*(Mr[i]-a)+(Mi[i]-b)*(Mi[i]-b));
        if ( distance <= min) {
                min = distance;
                kkk = i;
        }
    }
    err[0] = Mr[kkk] - a; 
    err[1] = Mi[kkk] - b;
}




/* Updating of coefficients */
void updatecoeff(double *hr,double *hi, int Ntap,
                 double *xr,double *xi, int Ndim,
                 double outr, double outi, double mu, 
                 double *Mr,double *Mi, int N)
{
    double *err;
    int i;
    err = mxMalloc(2*sizeof(double));
    errorfun( outr, outi, Mr, Mi, N, err);
    for(i=0;i<Ntap;i++) {
        *(hr+i) += mu*( err[0]*xr[i] + err[1]*xi[i]);
        *(hi+i) += mu*( err[1]*xr[i] - err[0]*xi[i]);
        *(hr+i+Ntap) += mu*( err[0]*xr[i+Ndim] + err[1]*xi[i+Ndim]);
        *(hi+i+Ntap) += mu*( err[1]*xr[i+Ndim] - err[0]*xi[i+Ndim]);
    }
    mxFree(err);
}



/* Core LMS adaptive filter function */
void lmsfilter(double *xr, double *xi, int Ndim, 
               double *h1r,double *h1i,double *h2r,double *h2i, int Ntap,
               double mu,double *Mr,double *Mi,
               double *yr,double *yi,double *mse,double *det,
               double sps,
               bool dontskip, int N)
{
    int i,j;
    int halftap = (Ntap-1)/2;
    int Ydim = Ndim-Ntap+1;
    double sum = 0, real, imag;
    for(i=0;i<Ydim;i++)
    {
        *(yr+i) = vmac(xr+i,h1r,Ntap) - vmac(xi+i,h1i,Ntap)
               +vmac(xr+i+Ndim,h1r+Ntap,Ntap)-vmac(xi+i+Ndim,h1i+Ntap,Ntap); 
        *(yi+i) = vmac(xi+i,h1r,Ntap) + vmac(xr+i,h1i,Ntap) 
               +vmac(xi+i+Ndim,h1r+Ntap,Ntap)+vmac(xr+i+Ndim,h1i+Ntap,Ntap); 
        *(yr+i+Ydim) = vmac(xr+i,h2r,Ntap) - vmac(xi+i,h2i,Ntap) 
               +vmac(xr+i+Ndim,h2r+Ntap,Ntap)-vmac(xi+i+Ndim,h2i+Ntap,Ntap); 
        *(yi+i+Ydim) = vmac(xi+i,h2r,Ntap) + vmac(xr+i,h2i,Ntap) 
               +vmac(xi+i+Ndim,h2r+Ntap,Ntap)+vmac(xr+i+Ndim,h2i+Ntap,Ntap); 
        
        /* Updating filter coefficients */
        if (  dontskip || ( i % (int)sps == 0) )
        {
            updatecoeff(h1r,h1i,Ntap,xr+i,xi+i,Ndim,*(yr+i),     *(yi+i),     mu,Mr,Mi,N);
            updatecoeff(h2r,h2i,Ntap,xr+i,xi+i,Ndim,*(yr+i+Ydim),*(yi+i+Ydim),mu,Mr,Mi,N);
        }

        /* Calculate MSE and determinate of H */
        sum += 0;
        *(mse+i) = sum / (i+1);
        
        real = (*(h1r+halftap)* *(h2r+Ntap+halftap)-*(h1i+halftap)* *(h2i+Ntap+halftap))
			-(*(h1r+Ntap+halftap)* *(h2r+halftap)-*(h1i+Ntap+halftap)* *(h2i+halftap));
        imag = (*(h1r+halftap)* *(h2i+Ntap+halftap)+*(h1i+halftap)* *(h2r+Ntap+halftap))
			-(*(h1r+Ntap+halftap)* *(h2i+halftap)+*(h1i+Ntap+halftap)* *(h2r+halftap));
        
        *(det+i) = Q_sqrt(real*real+imag*imag);        
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,  const mxArray *prhs[])
{
    double *xr		= mxGetPr(prhs[0]);         // Input field vector, real part
    double *xi		= mxGetPi(prhs[0]);         // Input field vector, imag part
    double *h1r		= mxGetPr(prhs[1]);         // 1st output filter, real part
    double *h1i		= mxGetPi(prhs[1]);         // 1st output filter, imag part
    double *h2r		= mxGetPr(prhs[2]);         // 2nd output filter, real part
    double *h2i		= mxGetPi(prhs[2]);         // 2nd output filter, imag part
    double Ntap		= mxGetScalar(prhs[3]);     // Length of filters, real scalar
    double mu		= mxGetScalar(prhs[4]);     // Algo constant, real scalar
    double *Mr		= mxGetPr(prhs[5]);         // Constellation points
	double *Mi		= mxGetPi(prhs[5]);
    double sps		= mxGetScalar(prhs[6]);     // Samples x symbol, 1 or 2
    
    int Mdim		= mxGetM(prhs[0]);			// Length of input field vector
    int Npol		= mxGetN(prhs[0]);			// Should be 2
    int Mfilter1	= mxGetM(prhs[1]);			// Should be == Ntap
    int Nfilter1	= mxGetN(prhs[1]);			// Should be 2
    int Mfilter2	= mxGetM(prhs[2]);			// Should be == Ntap
    int Nfilter2	= mxGetN(prhs[2]);			// Should be 2
    int N1			= mxGetM(prhs[5]);               
    int N2			= mxGetN(prhs[5]);
    int Nconstell	= N1 * N2;					// # of constellation points
    
    double *yr,*yi,*MSE,*DET;
    
    bool dontskip;    

    if ((int)Ntap % 2 == 0)
        mexErrMsgTxt("Ntaps should be an ODD INTEGER.");
    if (!(Ntap == Mfilter1)||!(Ntap == Mfilter2))
        mexErrMsgTxt("Filter length should be equal.");
    if (!(Npol == 2)||!(Nfilter1 == 2)||!(Nfilter2 == 2))
        mexErrMsgTxt("Two polarizations needed.");
        
    /* Chechking the value of samples per symbol    */
    if ((int)sps == 1)
        dontskip = 1;
    else
        dontskip = 0;

    /* Chech if the vector that are supposed to be complex are really complex */
    /* If X, H1 or H2 are purely real numbers (which is often the case for H1 */
    /* and H2) then matlab does not allocates the memory for the imaginary    */
    /* part and the program will segfault if it tries to access it            */
    
    if (xi==NULL) {
        xi = mxCalloc(Mdim*Npol, sizeof(double));
        mxSetPi( (mxArray *)prhs[0], xi);
    }
    if (h1i==NULL) {
        h1i = mxCalloc(Ntap*Npol,sizeof(double));
        mxSetPi( (mxArray *)prhs[1], h1i);
    }
    if (h2i==NULL) { 
        h2i = mxCalloc(Ntap*Npol,sizeof(double));
        mxSetPi( (mxArray *)prhs[2], h2i);
    }
        
	/* Next line allocates memory for a (Mdim-Ntap+1) x Npol complex vector */
    plhs[0] = mxCreateDoubleMatrix(Mdim-Ntap+1,Npol,mxCOMPLEX);
	yr = mxGetPr(plhs[0]);
    yi = mxGetPi(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(Mdim-Ntap+1,1,mxREAL);
    MSE = mxGetPr(plhs[1]);

    plhs[2] = mxCreateDoubleMatrix(Mdim-Ntap+1,1,mxREAL);
    DET = mxGetPr(plhs[2]);
    
    lmsfilter(xr,xi,Mdim,h1r,h1i,h2r,h2i,Ntap,mu,Mr,Mi,yr,yi,MSE,DET,sps,dontskip,Nconstell);

    return;
}
