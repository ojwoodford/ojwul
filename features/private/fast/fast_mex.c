#include "mex.h"
#include "fast.h"
#include <stdint.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    xy* corners;
    int numcorners, xsize, ysize, a, threshold, type = 9;
    const uint8_t* data;
    uint16_t* out;
    
    /* Parse the arguments */
    if (nlhs != 1 || nrhs < 2 || nrhs > 3)
        mexErrMsgTxt("Unexpected arguments");
    if (mxGetClassID(prhs[0]) != mxUINT8_CLASS || mxIsComplex(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("data must be a real uint8 matrix");
    data = (const uint8_t*)mxGetData(prhs[0]);
    xsize = mxGetM(prhs[0]);
    ysize = mxGetN(prhs[0]);
    if (mxGetNumberOfElements(prhs[1]) != 1 || mxIsComplex(prhs[1]))
        mexErrMsgTxt("threshold not a valid value");
    threshold = (int)mxGetScalar(prhs[1]);
    if (nrhs > 2) {
        if (mxGetNumberOfElements(prhs[2]) != 1 || mxIsComplex(prhs[2]))
            mexErrMsgTxt("type not a valid value");
        type = (int)mxGetScalar(prhs[2]);
    }
    
    /* Do the corner detection */
    switch (type) {
        case 9:
            corners = fast9_detect_nonmax(data, xsize, ysize, xsize, threshold, &numcorners);
            break;
        case 10:
            corners = fast10_detect_nonmax(data, xsize, ysize, xsize, threshold, &numcorners);
            break;
        case 11:
            corners = fast11_detect_nonmax(data, xsize, ysize, xsize, threshold, &numcorners);
            break;
        case 12:
            corners = fast12_detect_nonmax(data, xsize, ysize, xsize, threshold, &numcorners);
            break;
        default:
            mexErrMsgTxt("type not a valid value");
            break;
    }
    
    /* Construct the output array */
    plhs[0] = mxCreateNumericMatrix(2, numcorners, mxUINT16_CLASS, mxREAL);
    out = (uint16_t *)mxGetData(plhs[0]);
    for (a = 0; a < numcorners; ++a) {
        out[a*2+0] = (uint16_t)corners[a].y+1;
        out[a*2+1] = (uint16_t)corners[a].x+1;
    }
    free(corners);
}
