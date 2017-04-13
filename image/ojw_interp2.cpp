// B = ojw_interp2(A, X, Y[, method[, oobv]])

// Written by ojw 20/9/06

#include <mex.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "interp2_methods.hpp"

// Define types
#include <stdint.h>

#if _MATLAB_ < 805 // R2015a
extern "C" mxArray *mxCreateUninitNumericArray(mwSize ndim, const size_t *dims, mxClassID classid, mxComplexity ComplexFlag);
#endif

template<class Method, class U, class V> static inline void wrapper_func4(U *B, U *G, Method &im, const V *X, const V *Y, const int num_points)
{
	// For each of the interpolation points
    int i;
    if (G == NULL) {
#pragma omp parallel for if (num_points > 1000) num_threads(2) default(shared) private(i)
        for (i = 0; i < num_points; ++i)
            im.lookup(&B[i], Y[i]-1.0, X[i]-1.0, num_points); // Do the interpolation
    } else {
#pragma omp parallel for if (num_points > 1000) num_threads(2) default(shared) private(i)
        for (i = 0; i < num_points; ++i)
            im.lookup_grad(&B[i], &G[i*2], Y[i]-1.0, X[i]-1.0, num_points, num_points*2); // Do the interpolation
    } 
	return;
}

template<class T, class U, class V> static inline void wrapper_func3(U *B, U *G, const T *A, const V *X, const V *Y, const int num_points, const int w, const int h, const int col, const U oobv, const char meth)
{
    // Call the fourth wrapper function according to the interpolation type
    switch (meth) {
        case 'l':
        {
            IM_LIN<T,U,V> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points);
            break;
        }
        case 's':
        {
            IM_SHIFT<T,U,V> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points);
            break;
        }
        case 'n':
        {
            if (G != NULL)
                mexErrMsgTxt("Gradient computation not supported for nearest interpolation method");
            IM_NEAR<T,U,V> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points);
            break;
        }
        case 'c':
        {
            if (G != NULL)
                mexErrMsgTxt("Gradient computation not supported for cubic interpolation method");
            IM_CUB<T,U,V> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points);
            break;
        }
        default:
            mexErrMsgTxt("Unsupported interpolation method");
            break;
    }
}

template<class T, class V> static inline void wrapper_func2(void *B, void *G, const T *A, const V *X, const V *Y, const int num_points, const int w, const int h, const int col, const double oobv, const char meth, const mxClassID out_class)
{
	// Call the third wrapper function according to the output type
	switch (out_class) {
		case mxDOUBLE_CLASS:
            wrapper_func3((double *)B, (double *)G, A, X, Y, num_points, w, h, col, (const double)oobv, meth);
			break;
		case mxSINGLE_CLASS:
            wrapper_func3((float *)B, (float *)G, A, X, Y, num_points, w, h, col, (const float)oobv, meth);
			break;
		case mxUINT8_CLASS:
            wrapper_func3((uint8_t *)B, (uint8_t *)G, A, X, Y, num_points, w, h, col, (const uint8_t)oobv, meth);
			break;
		case mxINT16_CLASS:
            wrapper_func3((int16_t *)B, (int16_t *)G, A, X, Y, num_points, w, h, col, (const int16_t)oobv, meth);
			break;
		case mxUINT16_CLASS:
            wrapper_func3((uint16_t *)B, (uint16_t *)G, A, X, Y, num_points, w, h, col, (const uint16_t)oobv, meth);
			break;
		case mxLOGICAL_CLASS:
            wrapper_func3((mxLogical *)B, (mxLogical *)G, A, X, Y, num_points, w, h, col, (const mxLogical)oobv, meth);
			break;
		default:
			mexErrMsgTxt("Unsupported output type");
			break;
	}
	return;
}

template<class T> static inline void wrapper_func(void *B, void *G, const T *A, const mxArray *prhs[], const int num_points, const int w, const int h, const int col, const double oobv, const char meth, const mxClassID out_class, const mxClassID in_class)
{
    // Get pointers to the coordinate arrays
    const void *X = mxGetData(prhs[1]);
    const void *Y = mxGetData(prhs[2]);
    
	// Call the second wrapper function according to the coordinate type
	switch (in_class) {
		case mxDOUBLE_CLASS:
			wrapper_func2(B, G, A, (const double *)X, (const double *)Y, num_points, w, h, col, oobv, meth, out_class);
			break;
		case mxSINGLE_CLASS:
			wrapper_func2(B, G, A, (const float *)X, (const float *)Y, num_points, w, h, col, oobv, meth, out_class);
			break;
		default:
			mexErrMsgTxt("X and Y must be floating point arrays");
			break;
	}
	return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs < 3 || nrhs > 5)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs < 1 || nlhs > 2)
		mexErrMsgTxt("Unexpected number of output arguments.");

	// Check argument types are valid
	mxClassID in_class = mxGetClassID(prhs[1]);
	if (mxGetClassID(prhs[2]) != in_class)
		mexErrMsgTxt("X and Y must be arrays of the same type");
	for (int i = 0; i < nrhs; ++i) {
        if (mxIsComplex(prhs[i]))
			mexErrMsgTxt("Inputs cannot be complex.");
	}

	// Get and check array dimensions
	int num_points = mxGetNumberOfElements(prhs[1]);
	if (num_points != mxGetNumberOfElements(prhs[2]))
		mexErrMsgTxt("X and Y must have the same dimensions");
	int ndims = mxGetNumberOfDimensions(prhs[0]);
	if (ndims > 20)
		mexErrMsgTxt("A has an unsupported number of dimensions");
	size_t out_dims[21];
    out_dims[0] = 2;
	out_dims[1] = mxGetM(prhs[1]);
	out_dims[2] = mxGetN(prhs[1]);
	const mwSize *dims = mxGetDimensions(prhs[0]);
	out_dims[3] = 1;
	int nchannels = 1;
	for (int i = 2; i < ndims; ++i) {
		out_dims[i+1] = dims[i];
		nchannels *= dims[i];
	}
	
	// Get the out of bounds value (oobv) and set the output class to the same class as the oobv.
	double oobv;
	mxClassID out_class;
	if (nrhs > 4) {
		// Get the value for oobv
		if (mxGetNumberOfElements(prhs[4]) != 1)
			mexErrMsgTxt("oobv must be a scalar.");
		oobv = mxGetScalar(prhs[4]);
		out_class = mxGetClassID(prhs[4]);
	} else {
		// Use the default value for oobv
		oobv = mxGetNaN();
		out_class = in_class;
	}
	
	// Create the output arrays
	plhs[0] = mxCreateUninitNumericArray(ndims, &out_dims[1], out_class, mxREAL);
	void *B = mxGetData(plhs[0]);
    void *G = NULL;
    if (nlhs > 1) {
        plhs[1] = mxCreateUninitNumericArray(ndims+1, out_dims, out_class, mxREAL);
        G = mxGetData(plhs[1]);
    }

	// Get the interpolation method
    char buffer[10] = {'l'};
    int k = 0;
    if (nrhs > 3) {
        // Read in the method string
        if (mxGetString(prhs[3], buffer, sizeof(buffer)))
            mexErrMsgTxt("Unrecognised interpolation method");
        // Remove '*' from the start
        k += (buffer[k] == '*');
        // Remove 'bi' from the start
        k += 2 * ((buffer[k] == 'b') & (buffer[k+1] == 'i'));
    }
    
    // Get pointer to the input image
    const void *A = mxGetData(prhs[0]);
    
	// Call the first wrapper function according to the input image type
	switch (mxGetClassID(prhs[0])) {
		case mxDOUBLE_CLASS:
			wrapper_func(B, G, (const double *)A, prhs, num_points, dims[1], dims[0], nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxSINGLE_CLASS:
			wrapper_func(B, G, (const float *)A, prhs, num_points, dims[1], dims[0], nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxUINT8_CLASS:
			wrapper_func(B, G, (const uint8_t *)A, prhs, num_points, dims[1], dims[0], nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxINT16_CLASS:
			wrapper_func(B, G, (const int16_t *)A, prhs, num_points, dims[1], dims[0], nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxUINT16_CLASS:
			wrapper_func(B, G, (const uint16_t *)A, prhs, num_points, dims[1], dims[0], nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxLOGICAL_CLASS:
			wrapper_func(B, G, (const mxLogical *)A, prhs, num_points, dims[1], dims[0], nchannels, oobv, buffer[k], out_class, in_class);
			break;
		default:
			mexErrMsgTxt("A is of an unsupported type");
			break;
	}
	return;
}