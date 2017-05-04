#include <mex.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdint.h>
#include "interp2_methods.hpp"
#include "../utils/private/class_handle.hpp"
typedef IM_SHIFT<uint8_t, double, double> shiftIm_t;

#if _MATLAB_ < 805 // R2015a
extern "C" mxArray *mxCreateUninitNumericArray(mwSize ndim, const size_t *dims, mxClassID classid, mxComplexity ComplexFlag);
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs < 1 || nrhs > 3)
		mexErrMsgTxt("Unexpected number of input arguments.");
    
    // Destructor
    if (nrhs == 1) {
        // Destroy the C++ object
        destroyObject<shiftIm_t>(prhs[0]);
        if (nlhs != 0)
            mexWarnMsgTxt("Unexpected number of output arguments.");
        return;
    }
    
	if (nlhs != 1)
		mexErrMsgTxt("Unexpected number of output arguments.");
    
    // Constructor
    if (nrhs == 2) {
        // Get the image matrix
        if (mxGetClassID(prhs[0]) != mxUINT8_CLASS)
            mexErrMsgTxt("im must be a uint8 array.");
        const uint8_t *im = (const uint8_t *)mxGetData(prhs[0]);
        // Get the number of channels
        int ndims = mxGetNumberOfDimensions(prhs[0]);
        const int *dims = mxGetDimensions(prhs[0]);
        int nchannels = 1;
        for (int i = 2; i < ndims; ++i)
            nchannels *= dims[i];
        
        // Get the value for oobv
        if (mxGetNumberOfElements(prhs[1]) != 1 || !mxIsDouble(prhs[1]))
            mexErrMsgTxt("oobv must be a scalar double.");
        double oobv = mxGetScalar(prhs[1]);
        
        // Create a new instance of the interpolator
        shiftIm_t *instance = new shiftIm_t(im, oobv, dims[0], dims[1], nchannels);
        
        // Return a handle to the class
        plhs[0] = convertPtr2Mat<shiftIm_t>(instance);
        return;
    }
    
    // Get the class instance pointer from the first input
    shiftIm_t *instance = convertMat2Ptr<shiftIm_t>(prhs[0]);

	// Check argument types are valid
    int num_points = mxGetNumberOfElements(prhs[1]);
	if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || num_points != mxGetNumberOfElements(prhs[2]) || mxIsComplex(prhs[1]) || mxIsComplex(prhs[2]))
		mexErrMsgTxt("X and Y must be real double arrays of the same size");
    
    // Get pointers to the coordinate arrays
    const double *X = (const double *)mxGetData(prhs[1]);
    const double *Y = (const double *)mxGetData(prhs[2]);

	// Get output dimensions
	size_t out_dims[3];
	out_dims[0] = mxGetM(prhs[1]);
	out_dims[1] = mxGetN(prhs[1]);
    out_dims[2] = instance->Channels();
    
    // Create the output array
    plhs[0] = mxCreateUninitNumericArray(3, out_dims, mxDOUBLE_CLASS, mxREAL);
    double *out = (double *)mxGetData(plhs[0]);
    
    // For each of the interpolation points
    int i;
#pragma omp parallel for if (num_points > 1000) num_threads(2) default(shared) private(i)
	for (i = 0; i < num_points; ++i)
        instance->lookup(&out[i], Y[i]-1.0, X[i]-1.0, num_points); // Do the interpolation
	
	return;
}