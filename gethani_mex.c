#include "mex.h"
#include <math.h>

extern void cgethani(double x, double y, int kern, double value, double wght[2], double eps, double *bw);

// Gateway function gethani 
//  for function cgethani in gethani.c
//  to be called:
//  hakt = gethani(x,      # 1
//                 y,      # 2
//                 kern,   # 3
//                 value,  # 4
//                 wght,   # 5
//                 eps);   # 6
//

// MEX function gateway
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double x, y, value, eps, *wght, bw;
    int kern;
    
    // Check for proper number of arguments
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("MATLAB:gethani:invalidNumInputs", "GETHANI requires 6 input arguments");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:gethani:invalidNumOutputs", "GETHANI requires 1 output argument");
    }
    
    // Check if inputs are numeric
    if (!mxIsNumeric(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB:gethani:invalidInput", "X has to be a number");
    }
    if (!mxIsNumeric(prhs[1])) {
        mexErrMsgIdAndTxt("MATLAB:gethani:invalidInput", "Y has to be a number");
    }
    if (!mxIsNumeric(prhs[2])) {
        mexErrMsgIdAndTxt("MATLAB:gethani:invalidInput", "KERN has to be a number");
    }
    if (!mxIsNumeric(prhs[3])) {
        mexErrMsgIdAndTxt("MATLAB:gethani:invalidInput", "VALUE has to be a number");
    }
    if (!mxIsNumeric(prhs[4])) {
        mexErrMsgIdAndTxt("MATLAB:gethani:invalidInput", "WGHT has to be a vector");
    }
    if (!mxIsNumeric(prhs[5])) {
        mexErrMsgIdAndTxt("MATLAB:gethani:invalidInput", "EPS has to be a number");
    }
    
    // Get scalar inputs
    x = mxGetScalar(prhs[0]);
    y = mxGetScalar(prhs[1]);
    kern = (int)mxGetScalar(prhs[2]);
    value = mxGetScalar(prhs[3]);
    eps = mxGetScalar(prhs[5]);
    
    // Get wght vector
    if (mxGetNumberOfElements(prhs[4]) != 2) {
        mexErrMsgIdAndTxt("MATLAB:gethani:invalidInput", "WGHT has to be a vector of length 2");
    }
    wght = mxGetPr(prhs[4]);
    
    // Call the gethani function
    cgethani(x, y, kern, value, wght, eps, &bw);
    
    // Create output
    plhs[0] = mxCreateDoubleScalar(bw);
}