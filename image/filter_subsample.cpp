// B = filter_subsample(A, F)

// $Id: filter_subsample.cpp,v 1.4 2007/03/15 14:08:00 ojw Exp $

#include <string.h>
#include <mex.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_procs()  1
#define omp_get_thread_num() 0
#endif

#define MAX_NCHANNELS 10         // Allocate max number of channel pointers on the stack - alternative: use 'new'

#if _MATLAB_ < 805 // R2015a
extern "C" mxArray *mxCreateUninitNumericArray(mwSize ndim, const size_t *dims, mxClassID classid, mxComplexity ComplexFlag);
#endif

template<class T, class U> static inline void filter_subsample(const mxArray *prhs[], mxArray *plhs[], mxClassID in_class);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs != 2)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Unexpected number of output arguments.");

	// Check argument types are valid
	mxClassID in_class = mxGetClassID(prhs[0]);
	if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS && (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS || in_class != mxDOUBLE_CLASS))
		mexErrMsgTxt("F must be singles, or doubles iff A is doubles");
	if (mxGetM(prhs[1]) != 1 && mxGetN(prhs[1]) != 1)
		mexErrMsgTxt("F must be a vector");
	for (int i = 0; i < nrhs; i++) {
        if (mxIsComplex(prhs[i]))
			mexErrMsgTxt("Inputs cannot be complex.");
	}

	// Call the wrapper function according to the input type
	switch (in_class) {
		case mxLOGICAL_CLASS:
			filter_subsample<mxLogical, float>(prhs, plhs, in_class);
			break;
		case mxUINT8_CLASS:
			filter_subsample<unsigned char, float>(prhs, plhs, in_class);
			break;
		case mxUINT16_CLASS:
			filter_subsample<unsigned short, float>(prhs, plhs, in_class);
			break;
		case mxSINGLE_CLASS:
			filter_subsample<float, float>(prhs, plhs, in_class);
			break;
		case mxDOUBLE_CLASS:
			filter_subsample<double, double>(prhs, plhs, in_class);
			break;
		default:
			mexErrMsgTxt("A is of an unsupported type");
			break;
	}
	return;
}

static inline int sp(int x, int w) { w--; return abs(((x + (w * 9)) % (2 * w)) - w); } // Macro for symmetric padding

// Function for correct rounding
// Add these to use numeric_limits class
#include <limits>
using namespace std;
template<class U, class T> static inline U saturate_cast(T val)
{
	if (numeric_limits<U>::is_integer && !numeric_limits<T>::is_integer) {
		if (numeric_limits<U>::is_signed)
			return val > 0 ? (val > (T)numeric_limits<U>::max() ? numeric_limits<U>::max() : static_cast<U>(val + 0.5)) : (val < (T)numeric_limits<U>::min() ? numeric_limits<U>::min() : static_cast<U>(val - 0.5));
		else
			return val > 0 ? (val > (T)numeric_limits<U>::max() ? numeric_limits<U>::max() : static_cast<U>(val + 0.5)) : 0;
	}
	return static_cast<U>(val);
}

template<class T, class U> static inline void filter_subsample(const mxArray *prhs[], mxArray *plhs[], mxClassID in_class)
{
	// Get input data and sizes
	const mwSize *in_dims = mxGetDimensions(prhs[0]);
	size_t out_dims[3];
	out_dims[0] = (in_dims[0] + 1) / 2;
	out_dims[1] = (in_dims[1] + 1) / 2;
	out_dims[2] = mxGetNumberOfElements(prhs[0]) / (in_dims[0] * in_dims[1]);
	if (out_dims[2] > MAX_NCHANNELS)
		mexErrMsgTxt("A has more than the maximum MAX_NCHANNELS channels");
	int filter_length = mxGetNumberOfElements(prhs[1]);
	int col_pitch = out_dims[0] + filter_length;
	int co = filter_length / 2;
	U *F = (U *)mxGetData(prhs[1]);

	// Create the output and temporary arrays
	plhs[0] = mxCreateUninitNumericArray(3, out_dims, in_class, mxREAL);
	U *tmp_buf = (U *)mxCalloc((in_dims[1]+filter_length)*col_pitch*out_dims[2], sizeof(U));

	// Initialize pointers
	const T *A_ = (const T *)mxGetData(prhs[0]);
    U *D_ = tmp_buf + co + 1;
    int Aoff[MAX_NCHANNELS];
    int Doff[MAX_NCHANNELS];
    Aoff[0] = 0;
    Doff[0] = 0;
	for (int chan = 1; chan < out_dims[2]; ++chan) {
		Aoff[chan] = Aoff[chan-1] + in_dims[0] * in_dims[1];
		Doff[chan] = Doff[chan-1] + (in_dims[1] + filter_length) * col_pitch;
	}

	// Filter and subsample columns
	// For each column...
#pragma omp parallel for num_threads(omp_get_num_procs())
	for (int c = 0; c < in_dims[1]; ++c) {
        // Create pointers to new column
		int start = 1;
		int r = -co;
		int r2 = 0;
        const T *A = A_ + in_dims[0] * c;
        U *D = D_ + col_pitch * c;
		// For each row...
		// Boundary pixels
		for ( ; r < 0; ++r, r2 += start) {
			start = 1 - start;
			// For each channel...
			int index = sp(r,in_dims[0]); // Account for out of bound pixels
			for (int chan = 0; chan < out_dims[2]; ++chan) {
				U val = (U)A[index+Aoff[chan]];
				// For each filter entry...
				for (int f = start, i = 0; f < filter_length; f += 2, ++i)
					D[(Doff[chan]+r2)-i] +=  F[f] * val;
			}
		}
		// Centre pixels
		for ( ; r < in_dims[0]; ++r, r2 += start) {
			start = 1 - start;
			// For each channel...
			for (int chan = 0; chan < out_dims[2]; ++chan) {
				U val = (U)A[r+Aoff[chan]];
				// For each filter entry...
				for (int f = start, i = 0; f < filter_length; f += 2, ++i)
					D[(Doff[chan]+r2)-i] +=  F[f] * val;
			}
		}
		// Boundary pixels
		for ( ; r < in_dims[0]+filter_length-co-1; ++r, r2 += start) {
			start = 1 - start;
			// For each channel...
			int index = sp(r,in_dims[0]); // Account for out of bound pixels
			for (int chan = 0; chan < out_dims[2]; ++chan) {
				U val = (U)A[index+Aoff[chan]];
				// For each filter entry...
				for (int f = start, i = 0; f < filter_length; f += 2, ++i)
					D[(Doff[chan]+r2)-i] +=  F[f] * val;
			}
		}
	}

	// Initialize pointers
    T *O_ = (T *)mxGetData(plhs[0]);
    int Boff[MAX_NCHANNELS];
    Boff[0] = 0;
    for (int chan = 1; chan < out_dims[2]; ++chan) {
		Aoff[chan] = Aoff[chan-1] + out_dims[0] * out_dims[1];
		Boff[chan] = Boff[chan-1] + out_dims[1] + filter_length;
    }
    int row_pitch = (out_dims[1] + filter_length) * out_dims[2];
	U *tmp_buf2 = (U *)mxMalloc(row_pitch*omp_get_num_procs()*sizeof(U));
	U *B_ = tmp_buf2 + co + 1;

	// Filter and subsample rows
	// For each row...
#pragma omp parallel for num_threads(omp_get_num_procs())
	for (int r = 0; r < out_dims[0]; ++r) {
        // Create pointers to new row
		int start = 1;
		int c2 = 0;
		int c = -co;
        U *D = D_ + r;
        U *B = B_ + row_pitch * omp_get_thread_num();
		// Clear the row buffer
        memset(tmp_buf2+row_pitch*omp_get_thread_num(), 0, row_pitch*sizeof(U));
		// For each column...
		// Boundary pixels
		for ( ; c < 0; ++c, c2 += start) {
			start = 1 - start;
			// For each channel...
			int index = sp(c,in_dims[1]); // Account for out of bound pixels
			for (int chan = 0; chan < out_dims[2]; ++chan) {
				U val = D[index*col_pitch+Doff[chan]];
				// For each filter entry...
				for (int f = start, i = 0; f < filter_length; f += 2, ++i)
                    B[Boff[chan]+(c2-i)] +=  F[f] * val;
			}
		}
		// Centre pixels
		for ( ; c < in_dims[1]; ++c, c2 += start) {
			start = 1 - start;
			// For each channel...
			for (int chan = 0; chan < out_dims[2]; ++chan) {
				U val = D[c*col_pitch+Doff[chan]];
				// For each filter entry...
				for (int f = start, i = 0; f < filter_length; f += 2, ++i)
                    B[Boff[chan]+(c2-i)] +=  F[f] * val;
			}
		}
		// Boundary pixels
		for ( ; c < in_dims[1]+filter_length-co-1; ++c, c2 += start) {
			start = 1 - start;
			// For each channel...
			int index = sp(c,in_dims[1]); // Account for out of bound pixels
			for (int chan = 0; chan < out_dims[2]; ++chan) {
				U val = D[index*col_pitch+Doff[chan]];
				// For each filter entry...
				for (int f = start, i = 0; f < filter_length; f += 2, ++i)
                    B[Boff[chan]+(c2-i)] +=  F[f] * val;
			}
		}
        // Copy to the output buffer
        T *O = O_ + r;
        for (c = 0; c < out_dims[1]; ++c) {
            for (int chan = 0; chan < out_dims[2]; ++chan)
                O[c*out_dims[0]+Aoff[chan]] = saturate_cast<T, U>(B[c+Boff[chan]]);
        }
	}

    // Free the temp buffers
	mxFree(tmp_buf2);
	mxFree(tmp_buf);
	return;
}
