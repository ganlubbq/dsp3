#include <math.h>
#include "mex.h"

/* The Polytechnic University	version: 2011/09/07 office pc	*/
/* The Polytechnic University	version: 2012/11/14 office pc	*/

/* Multiply and Accumulate a vector */
double vmac(double *a, double *b, int n)
{
    double z = 0;
    int i=0;
    for(i=0;i<n;i++)
        z += a[i] * b[i];
    return z;
}



void errorfun(double a,double b, double *Mr, double *Mi, int N, double *err)
{
    double min, distance;
    int kkk, i;
    min = 10000000;  
    kkk = 0;
    for ( i=0;i<N;i++ ) {
        distance = sqrt((Mr[i]-a)*(Mr[i]-a)+(Mi[i]-b)*(Mi[i]-b));
        if (distance <= min) {
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
    errorfun(outr, outi, Mr, Mi, N, err);
    for(i=0;i<Ntap;i++) {
        *(hr+i) += mu * (err[0] * xr[i] + err[1] * xi[i]);
        *(hi+i) += mu * (err[1] * xr[i] - err[0] * xi[i]);
    }
    mxFree(err);
}



/* Core LMS adaptive filter function */
void lmsfilter(double *xr, double *xi, int Ndim, 
               double *hr,double *hi,int Ntap,
               double mu,double *Mr,double *Mi,
               double *yr,double *yi,double *mse,double *det,
               double sps,
               bool dontskip, int N)
{
    int i;
    int halftap = (Ntap-1)/2;
    int Ydim = Ndim-Ntap+1;
    double sum=0,real=0,imag=0;
    for(i=0;i<Ydim;i++)
    {
        *(yr+i) = vmac(xr+i,hr,Ntap) - vmac(xi+i,hi,Ntap); 
        *(yi+i) = vmac(xi+i,hr,Ntap) + vmac(xr+i,hi,Ntap); 
        
        /* Updating filter coefficients */
        if (  dontskip || ( i % (int)sps == 0) ) {
            updatecoeff(hr,hi,Ntap,xr+i,xi+i,Ndim,*(yr+i),*(yi+i),mu,Mr,Mi,N);
        }

        /* Calculate MSE and determinate of H */
        sum += 0;
        *(mse+i) = sum / (i+1);
        
        real = *(hr+halftap);
        imag = *(hi+halftap);
        
        *(det+i) = sqrt(real*real+imag*imag);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,  const mxArray *prhs[])
{
	double *yr,*yi,*xr,*xi,*hr,*hi,*Mr,*Mi,*MSE,*DET;
    double sps,mu;
	int Ntap,Mdim,Npol,Mfilt,Nfilt,N1,N2,Nconstell;
    bool dontskip;

	if(nrhs != 6) 
		mexErrMsgTxt("6 input required.");
    if(nlhs > 3) 
		mexErrMsgTxt("Too many output arguments.");

    xr			= mxGetPr(prhs[0]);				// Input field vector, real part
    xi			= mxGetPi(prhs[0]);				// Input field vector, imag part
    hr			= mxGetPr(prhs[1]);				// 1st output filter, real part
    hi			= mxGetPi(prhs[1]);				// 1st output filter, imag part
    Ntap		= mxGetScalar(prhs[2]);			// Length of filters, real scalar
    mu			= mxGetScalar(prhs[3]);			// Algo constant, real scalar
    Mr			= mxGetPr(prhs[4]);				// Constellation points
	Mi			= mxGetPi(prhs[4]);
	sps			= mxGetScalar(prhs[5]);			// Samples x symbol, 1 or 2
    
    Mdim		= mxGetM(prhs[0]);				// Length of input field vector
    Npol		= mxGetN(prhs[0]);				// Should be 1
    Mfilt		= mxGetM(prhs[1]);				// Should be == Ntap
    Nfilt		= mxGetN(prhs[1]);				// Should be 1
    N1			= mxGetM(prhs[4]);               
    N2			= mxGetN(prhs[4]);
    Nconstell	= N1 * N2;						// Length of constellation points

    if (!(Ntap % 2))
        mexErrMsgTxt("Ntaps should be an ODD INTEGER.");
    if (!(Ntap == Mfilt))
        mexErrMsgTxt("Filter length should be equal.");
    if (!(Npol == 1)||!(Nfilt == 1))
        mexErrMsgTxt("Single polarization is required.");
        
    /* Chechking the value of samples per symbol */
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
    if (hi==NULL) {
        hi = mxCalloc(Ntap*Npol,sizeof(double));
        mxSetPi( (mxArray *)prhs[2], hi);
    }
        
	/* Next line allocates memory for a (Mdim-Ntap+1) x Npol complex vector */
    plhs[0] = mxCreateDoubleMatrix(Mdim-Ntap+1,Npol,mxCOMPLEX);
	yr = mxGetPr(plhs[0]);
    yi = mxGetPi(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(Mdim-Ntap+1,1,mxREAL);
    MSE = mxGetPr(plhs[1]);

    plhs[2] = mxCreateDoubleMatrix(Mdim-Ntap+1,1,mxREAL);
    DET = mxGetPr(plhs[2]);
    
    lmsfilter(xr,xi,Mdim,hr,hi,Ntap,mu,Mr,Mi,yr,yi,MSE,DET,sps,dontskip,Nconstell);

    return;
}
