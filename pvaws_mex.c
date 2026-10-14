#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

// Forward declaration of the pvaws function from aws.c
extern void pvaws(
    double *y,     // 1
    int *pos,      // 2
    int nv,        // 3
    int nvd,       // 4
    int n1,        // 5
    int n2,        // 6
    int n3,        // 7
    double hakt,   // 8
    double lambda, // 9
    double *theta, // 10
    double *bi,    // 11
    double *bin,   // 12
    double *thnew, // 13
    double *invcov,// 14
    double spmin,  // 16
    double *lwght, // 17
    double *wght,  // 18
    double *swjy,  // 19
    int np1,       // 20
    int np2,       // 21
    int np3        // 22
);

//     Gateway function pvaws
//     for C function pvaws in aws.c
//     
//      [bi, theta] = pvaws(modelCoeff,              1
//                          position,                2
//                          nvec,                    3
//                          nvec * (nvec + 1) / 2,   4
//                          n1,                      5
//                          n2,                      6
//                          n3,                      7
//                          hakt,                    8
//                          lambda0,                 9
//                          theta,                  10
//                          bi,                     11
//                          invCov,                 12
//                          spmin,                  13 
//                          wghts,                  14
//                          dlw,                    15
//                          np1,                    16
//                          np2,                    17
//                          np3);                   18


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    // local integer variables

    mwSize NV_IN;
    mwSize NVD_IN;
    mwSize N1_IN;
    mwSize N2_IN;
    mwSize N3_IN;
    double SPMIN_IN;
    mwSize DLW_IN;
    mwSize NP1_IN;
    mwSize NP2_IN;
    mwSize NP3_IN;

    // local real variables (non-portable syntax)

    double HAKT_IN;
    double LAMBDA_IN;

    // local input pointers

    double *Y_IN_PR;
    int *POS_IN_PR;
    double *THETA_IN_PR;
    double *BI_IN_PR;
    double *BI_OUT_PR;
    double *THNEW_OUT_PR;
    double *INVCOV_IN_PR;
    double *LWGHT_LOCAL_PR;
    double *WGHT_IN_PR;
    double *SWJY_LOCAL_PR;

    // Check for proper number of arguments

    if (nrhs != 19)
        mexErrMsgTxt("pvaws requires 19 input arguments");
    if (nlhs != 2)
        mexErrMsgTxt("pvaws requires 2 output arguments");

    // get all the prhs

    Y_IN_PR = mxGetPr(prhs[0]);
    POS_IN_PR = (int*)mxGetData(prhs[1]);
    // Get scalar parameters with type checking
    if (!mxIsNumeric(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgTxt("NV_IN must be a scalar numeric value");
    }
    NV_IN = (mwSize)mxGetScalar(prhs[2]);

    if (!mxIsNumeric(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1) {
        mexErrMsgTxt("NVD_IN must be a scalar numeric value");
    }
    NVD_IN = (mwSize)mxGetScalar(prhs[3]);

    if (!mxIsNumeric(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1) {
        mexErrMsgTxt("N1_IN must be a scalar numeric value");
    }
    N1_IN = (mwSize)mxGetScalar(prhs[4]);

    if (!mxIsNumeric(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgTxt("N2_IN must be a scalar numeric value");
    }
    N2_IN = (mwSize)mxGetScalar(prhs[5]);

    if (!mxIsNumeric(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgTxt("N3_IN must be a scalar numeric value");
    }
    N3_IN = (mwSize)mxGetScalar(prhs[6]);
    HAKT_IN = mxGetScalar(prhs[7]);
    LAMBDA_IN = mxGetScalar(prhs[8]);
    THETA_IN_PR = mxGetPr(prhs[9]);
    BI_IN_PR = mxGetPr(prhs[10]);
    INVCOV_IN_PR = mxGetPr(prhs[11]);

    // prhs[12] = mc.cores (ncores) — ignored, C code runs single-threaded

    if (!mxIsNumeric(prhs[13]) || mxGetNumberOfElements(prhs[13]) != 1) {
        mexErrMsgTxt("SPMIN_IN must be a scalar numeric value");
    }
    SPMIN_IN = mxGetScalar(prhs[13]);

    WGHT_IN_PR = mxGetPr(prhs[14]);

    if (!mxIsNumeric(prhs[15]) || mxGetNumberOfElements(prhs[15]) != 1) {
        mexErrMsgTxt("DLW_IN must be a scalar numeric value");
    }
    DLW_IN = (mwSize)mxGetScalar(prhs[15]);

    if (!mxIsNumeric(prhs[16]) || mxGetNumberOfElements(prhs[16]) != 1) {
        mexErrMsgTxt("NP1_IN must be a scalar numeric value");
    }
    NP1_IN = (mwSize)mxGetScalar(prhs[16]);

    if (!mxIsNumeric(prhs[17]) || mxGetNumberOfElements(prhs[17]) != 1) {
        mexErrMsgTxt("NP2_IN must be a scalar numeric value");
    }
    NP2_IN = (mwSize)mxGetScalar(prhs[17]);

    if (!mxIsNumeric(prhs[18]) || mxGetNumberOfElements(prhs[18]) != 1) {
        mexErrMsgTxt("NP3_IN must be a scalar numeric value");
    }
    NP3_IN = (mwSize)mxGetScalar(prhs[18]);

    // create all the plhs

    mwSize nmask = mxGetNumberOfElements(prhs[10]);
    plhs[0] = mxCreateDoubleMatrix(1, nmask, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(NV_IN, nmask, mxREAL);

    if (plhs[0] == NULL || plhs[1] == NULL) {
        mexErrMsgTxt("Memory allocation failed for output arrays");
    }

    BI_OUT_PR = mxGetPr(plhs[0]);
    THNEW_OUT_PR = mxGetPr(plhs[1]);

    // create all arrays for internal calculations

    mxArray *LWGHT_LOCAL = mxCreateDoubleMatrix(1, DLW_IN, mxREAL);
    mxArray *SWJY_LOCAL = mxCreateDoubleMatrix(NV_IN, 1, mxREAL);

    if (LWGHT_LOCAL == NULL || SWJY_LOCAL == NULL) {
        mxDestroyArray(LWGHT_LOCAL);
        mxDestroyArray(SWJY_LOCAL);
        mexErrMsgTxt("Memory allocation failed for temporary arrays");
    }

    LWGHT_LOCAL_PR = mxGetPr(LWGHT_LOCAL);
    SWJY_LOCAL_PR = mxGetPr(SWJY_LOCAL);

     // Call the actual C function 'pvaws2' which performs the main computation for the pvaws algorithm.
     // This function processes the input data and computes the outputs 'BI_OUT_PR' and 'THNEW_OUT_PR'.

     pvaws(
         Y_IN_PR,           // 1
         POS_IN_PR,         // 2
         (int)NV_IN,        // 3
         (int)NVD_IN,       // 4
         (int)N1_IN,        // 5
         (int)N2_IN,        // 6
         (int)N3_IN,        // 7
         HAKT_IN,           // 8
         LAMBDA_IN,         // 9
         THETA_IN_PR,       // 10
         BI_IN_PR,          // 11
         BI_OUT_PR,         // 12
         THNEW_OUT_PR,      // 13
         INVCOV_IN_PR,      // 14
         SPMIN_IN,          // 15
         LWGHT_LOCAL_PR,    // 16
         WGHT_IN_PR,        // 17
         SWJY_LOCAL_PR,     // 18
         (int)NP1_IN,       // 19
         (int)NP2_IN,       // 20
         (int)NP3_IN        // 21
     );

    // Free allocated memory for temporary arrays used in internal calculations.
    // These arrays were created to store intermediate results and are no longer needed.
    mxDestroyArray(LWGHT_LOCAL);
    mxDestroyArray(SWJY_LOCAL);
}