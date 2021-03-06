// B = ojw_interp2(A, X, Y[, method[, oobv[, max_num_threads]]])

// Written by ojw 20/9/06

#include <mex.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "interp2_methods.hpp"
#include <vector>

// Define types
#include <stdint.h>

#if _MATLAB_ < 805 // R2015a
extern "C" mxArray *mxCreateUninitNumericArray(mwSize ndim, const size_t *dims, mxClassID classid, mxComplexity ComplexFlag);
#endif

template<class Method, class U, class V> static inline void wrapper_func4(U *B, U *G, Method &im, const V *X, const V *Y, const int num_points, const int num_threads)
{
	// For each of the interpolation points
#pragma omp parallel if (num_threads > 1) num_threads(num_threads) default(shared)
    {
        if (G == NULL) {
#pragma omp for schedule(dynamic, 128)
            for (int i = 0; i < num_points; ++i)
                im.lookup(&B[i], Y[i]-static_cast<V>(1.0), X[i]-static_cast<V>(1.0), num_points); // Do the interpolation
        } else {
#pragma omp for schedule(dynamic, 128)
            for (int i = 0; i < num_points; ++i)
                im.lookup_grad(&B[i], &G[i*2], Y[i]-static_cast<V>(1.0), X[i]-static_cast<V>(1.0), num_points, num_points*2); // Do the interpolation
        }
    }
    return;
}

template<class T, class U, class V> static inline void wrapper_func3(U *B, U *G, const T *A, const V *X, const V *Y, const int num_points, const int num_threads, const int w, const int h, const int col, const U oobv, const char meth)
{
    // Call the fourth wrapper function according to the interpolation type
    switch (meth) {
        case 'l':
        {
            IM_LIN<T,U,V> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points, num_threads);
            break;
        }
        case 's':
        {
            IM_SHIFT<T,U,V> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points, num_threads);
            break;
        }
        case 'n':
        {
            IM_NEAR<T,U,V> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points, num_threads);
            break;
        }
        case 'm':
        {
            IM_NTAP<T,U,V,3,magic> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points, num_threads);
            break;
        }
        case '4':
        {
            IM_NTAP<T,U,V,4,lanczos<4>> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points, num_threads);
            break;
        }
        case '6':
        {
            IM_NTAP<T,U,V,6,lanczos<6>> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points, num_threads);
            break;
        }
        case 'c':
        {
            if (G != NULL)
                mexErrMsgTxt("Gradient computation not supported for cubic interpolation method");
            IM_CUB<T,U,V> im(A, oobv, h, w, col);
			wrapper_func4(B, G, im, X, Y, num_points, num_threads);
            break;
        }
        default:
            mexErrMsgTxt("Unsupported interpolation method");
            break;
    }
}

template<class T, class V> static inline void wrapper_func2(void *B, void *G, const T *A, const V *X, const V *Y, const int num_points, const int num_threads, const int w, const int h, const int col, const double oobv, const char meth, const mxClassID out_class)
{
	// Call the third wrapper function according to the output type
	switch (out_class) {
		case mxDOUBLE_CLASS:
            wrapper_func3((double *)B, (double *)G, A, X, Y, num_points, num_threads, w, h, col, (const double)oobv, meth);
			break;
		case mxSINGLE_CLASS:
            wrapper_func3((float *)B, (float *)G, A, X, Y, num_points, num_threads, w, h, col, (const float)oobv, meth);
			break;
		case mxUINT8_CLASS:
            wrapper_func3((uint8_t *)B, (uint8_t *)G, A, X, Y, num_points, num_threads, w, h, col, (const uint8_t)oobv, meth);
			break;
		case mxINT8_CLASS:
            wrapper_func3((int8_t *)B, (int8_t *)G, A, X, Y, num_points, num_threads, w, h, col, (const int8_t)oobv, meth);
			break;
		case mxINT16_CLASS:
            wrapper_func3((int16_t *)B, (int16_t *)G, A, X, Y, num_points, num_threads, w, h, col, (const int16_t)oobv, meth);
			break;
		case mxUINT16_CLASS:
            wrapper_func3((uint16_t *)B, (uint16_t *)G, A, X, Y, num_points, num_threads, w, h, col, (const uint16_t)oobv, meth);
			break;
		case mxLOGICAL_CLASS:
            wrapper_func3((mxLogical *)B, (mxLogical *)G, A, X, Y, num_points, num_threads, w, h, col, (const mxLogical)oobv, meth);
			break;
		default:
			mexErrMsgTxt("Unsupported output type");
			break;
	}
	return;
}

template<class T> static inline void wrapper_func(void *B, void *G, const T *A, const mxArray *prhs[], const int num_points, const int num_threads, const int w, const int h, const int col, const double oobv, const char meth, const mxClassID out_class, const mxClassID in_class)
{
    // Get pointers to the coordinate arrays
    const void *X = mxGetData(prhs[1]);
    const void *Y = mxGetData(prhs[2]);
    
	// Call the second wrapper function according to the coordinate type
	switch (in_class) {
		case mxDOUBLE_CLASS:
			wrapper_func2(B, G, A, (const double *)X, (const double *)Y, num_points, num_threads, w, h, col, oobv, meth, out_class);
			break;
		case mxSINGLE_CLASS:
			wrapper_func2(B, G, A, (const float *)X, (const float *)Y, num_points, num_threads, w, h, col, oobv, meth, out_class);
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
	if (nrhs < 3 || nrhs > 6)
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
	int num_points = static_cast<int>(mxGetNumberOfElements(prhs[1]));
	if (num_points != static_cast<int>(mxGetNumberOfElements(prhs[2])))
		mexErrMsgTxt("X and Y must have the same dimensions");
	int ndims = static_cast<int>(mxGetNumberOfDimensions(prhs[0]));
	std::vector<size_t> out_dims(ndims+2);
    out_dims[0] = 2;
	out_dims[1] = mxGetM(prhs[1]);
	out_dims[2] = mxGetN(prhs[1]);
	const mwSize *dims = mxGetDimensions(prhs[0]);
	out_dims[3] = 1;
	int nchannels = 1;
	for (int i = 2; i < ndims; ++i) {
		out_dims[i+1] = dims[i];
		nchannels *= static_cast<int>(dims[i]);
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
        plhs[1] = mxCreateUninitNumericArray(ndims+1, &out_dims[0], out_class, mxREAL);
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
    
    // Get the maximum number of threads
    int num_threads = 1;
#ifdef _OPENMP
    num_threads += num_points >> 10;
    num_threads = std::min(num_threads, omp_get_num_procs());
    if (nrhs > 5) {
		if (mxGetNumberOfElements(prhs[4]) != 1)
			mexErrMsgTxt("max_num_threads must be a scalar.");
		num_threads = std::min(static_cast<int>(mxGetScalar(prhs[5])), num_threads);
    }
#endif
    
    // Get pointer to the input image
    const void *A = mxGetData(prhs[0]);
    
	// Call the first wrapper function according to the input image type
    const int width = static_cast<int>(dims[1]);
    const int height = static_cast<int>(dims[0]);
	switch (mxGetClassID(prhs[0])) {
		case mxDOUBLE_CLASS:
			wrapper_func(B, G, (const double *)A, prhs, num_points, num_threads, width, height, nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxSINGLE_CLASS:
			wrapper_func(B, G, (const float *)A, prhs, num_points, num_threads, width, height, nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxINT8_CLASS:
			wrapper_func(B, G, (const int8_t *)A, prhs, num_points, num_threads, width, height, nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxUINT8_CLASS:
			wrapper_func(B, G, (const uint8_t *)A, prhs, num_points, num_threads, width, height, nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxINT16_CLASS:
			wrapper_func(B, G, (const int16_t *)A, prhs, num_points, num_threads, width, height, nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxUINT16_CLASS:
			wrapper_func(B, G, (const uint16_t *)A, prhs, num_points, num_threads, width, height, nchannels, oobv, buffer[k], out_class, in_class);
			break;
		case mxLOGICAL_CLASS:
			wrapper_func(B, G, (const mxLogical *)A, prhs, num_points, num_threads, width, height, nchannels, oobv, buffer[k], out_class, in_class);
			break;
		default:
			mexErrMsgTxt("A is of an unsupported type");
			break;
	}
	return;
}
