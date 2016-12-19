#include <math.h>
#include "mex.h"

/* Fast exponential of a column real vector 		*/
/*							*/ 
/* fastexp(x) = exp(i*x) with x real column vector      */
/* 							*/
/* Author: Paolo Serena, 2009				*/
/* University of Parma, Italy				*/

/*
 *    This file is part of Optilux, the optical simulator toolbox.
 *    Copyright (C) 2009  Paolo Serena, <serena@tlc.unipr.it>
 *			 
 *    Optilux is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 * 
 *    Optilux is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* these 2 #define lines help make the later code more readable */
/* Input Arguments */
#define PARAMETER_IN    prhs[0]

/* Output Arguments */
#define RESULT_OUT  plhs[0]


void fastexp(double*yr, double*ypr, double*ypi,  int m) {
  while(m>0) {
    m--;
    ypr[m]=cos(yr[m]);
    ypi[m]=sin(yr[m]);
    
  }
}

void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray*prhs[] )
{ 
    double *ypr, *ypi; 
    double *t, *yr; 
    unsigned int m,n; 
    
    
    m = mxGetM(PARAMETER_IN); 
    n = mxGetN(PARAMETER_IN);
    
    /* Create a matrix for the return argument */ 
    RESULT_OUT = mxCreateDoubleMatrix(m, n, mxCOMPLEX); 
    
    /* Assign pointers to the various parameters */ 
    ypr = mxGetPr(RESULT_OUT);    
    ypi = mxGetPi(RESULT_OUT);    
    yr = mxGetPr(PARAMETER_IN);
        
    /* Do the actual computation*/
    fastexp(yr,ypr,ypi,m*n); 
    return;
}
