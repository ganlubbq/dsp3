/*=================================================================
 *This function is used to remove comma from a Matlab string, i.e.,
 *using space to replace comma (,)
 *
 *See Also: TekDPO72004.m
 *=============================================================*/
#include "mex.h"

void revord(char *input_buf, mwSize buflen, char *output_buf)
{
  mwSize i;

  if (buflen == 0) return;

  /* find comma, replace it */
  for(i=0;i<buflen-1;i++) 
      if (*(input_buf+i)==',')
          *(output_buf+i) = ' ';
      else
          *(output_buf+i) = *(input_buf+i);
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *input_buf, *output_buf;
    mwSize buflen;
    
    /* check for proper number of arguments */
    if(nrhs!=1) 
      mexErrMsgTxt("One input required.");
    else if(nlhs > 1) 
      mexErrMsgTxt("Too many output arguments.");

    /* input must be a string */
    if ( mxIsChar(prhs[0]) != 1)
      mexErrMsgTxt("Input must be a string.");

    /* input must be a row vector */
    if (mxGetM(prhs[0])!=1)
      mexErrMsgTxt("Input must be a row vector.");
    
    /* get the length of the input string */
    buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;

    /* allocate memory for output string */
    output_buf = mxCalloc(buflen, sizeof(char));

    /* copy the string data from prhs[0] into a C string input_ buf.    */
    input_buf = mxArrayToString(prhs[0]);
    
    if(input_buf == NULL) 
      mexErrMsgTxt("Could not convert input to string.");
    
    /* call the C subroutine */
    revord(input_buf, buflen, output_buf);

    /* set C-style string output_buf to MATLAB mexFunction output*/
    plhs[0] = mxCreateString(output_buf);
    mxFree(input_buf);
    return;
}

