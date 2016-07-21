// index = find_first(A, operator, value, start)

#include <mex.h>
#include <stdint.h>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

#define FUNC(NAME, OPERATOR) \
template<typename T> static inline uint32_t NAME(const T *A, const T val, const size_t numel, const uint32_t start) \
{ \
    uint32_t ind; \
    for (ind = start; ind < numel; ++ind) { if (A[ind] OPERATOR val) { return ind+1; }} \
    for (ind = 0; ind < start; ++ind) { if (A[ind] OPERATOR val) { return ind+1; }} \
    return 0; \
}

FUNC(neq, !=);
FUNC(eq, ==);
FUNC(lteq, <=);
FUNC(gteq, >=);
FUNC(lt, <);
FUNC(gt, >);

template<typename T> static inline void wrapper_func(const T *A, const char *op, const T *val, int vs, const uint32_t *start, int ss, const int M, const int N, uint32_t *out)
{   
    // Set the defaults
    T val_ = 0;
    if (val == NULL)
        val = &val_;
    uint32_t start_ = 0;
    if (start == NULL)
        start = &start_;
    
	// Call the second wrapper function according to the operator type
    int i;
	if (!strcmp(op, "~=")) {
#pragma omp parallel for if (N > 1 && N * M > 1e6) num_threads(2) default(shared) private(i)
        for (i = 0; i < N; ++i)
            out[i] = neq(&A[i*M], val[i&vs], M, start[i&ss]);
    } else if (!strcmp(op, "==")) {
#pragma omp parallel for if (N > 1 && N * M > 1e6) num_threads(2) default(shared) private(i)
        for (i = 0; i < N; ++i)
            out[i] = eq(&A[i*M], val[i&vs], M, start[i&ss]);
    } else if (!strcmp(op, "<=")) {
#pragma omp parallel for if (N > 1 && N * M > 1e6) num_threads(2) default(shared) private(i)
        for (i = 0; i < N; ++i)
            out[i] = lteq(&A[i*M], val[i&vs], M, start[i&ss]);
    } else if (!strcmp(op, ">=")) {
#pragma omp parallel for if (N > 1 && N * M > 1e6) num_threads(2) default(shared) private(i)
        for (i = 0; i < N; ++i)
            out[i] = gteq(&A[i*M], val[i&vs], M, start[i&ss]);
    } else if (!strcmp(op, "<")) {
#pragma omp parallel for if (N > 1 && N * M > 1e6) num_threads(2) default(shared) private(i)
        for (i = 0; i < N; ++i)
            out[i] = lt(&A[i*M], val[i&vs], M, start[i&ss]);
    } else if (!strcmp(op, ">")) {
#pragma omp parallel for if (N > 1 && N * M > 1e6) num_threads(2) default(shared) private(i)
        for (i = 0; i < N; ++i)
            out[i] = gt(&A[i*M], val[i&vs], M, start[i&ss]);
    } else
        mexErrMsgTxt("Comparison operator not recognized.");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs < 1 || nrhs > 4)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs > 1)
		mexErrMsgTxt("Unexpected number of output arguments.");
    
    // Set the default input arguments
    const void *val = NULL;
    const uint32_t *start = NULL;
    char op[3] = {'~', '=', 0};
    int vs = 0, ss = 0;

	// Check argument types are valid
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("A should be a real, numeric array");
	mxClassID in_class = mxGetClassID(prhs[0]);
    int M = mxGetM(prhs[0]);
    int N = mxGetN(prhs[0]);
    if (M == 1) {
        M = N;
        N = 1;
    }
    if (nrhs > 2) {
        if (mxGetString(prhs[1], op, 3))
            mexErrMsgTxt("operator is of unexpected type or size");
        vs = mxGetNumberOfElements(prhs[2]);
        if (mxGetClassID(prhs[2]) != in_class || mxIsComplex(prhs[2]) || (vs != 1 && vs != N))
            mexErrMsgTxt("value is of unexpected type or size");
        vs = vs == 1 ? 0 : -1;
        val = mxGetData(prhs[2]);
    }
    if (nrhs == 2 || nrhs == 4) {
        ss = mxGetNumberOfElements(prhs[nrhs-1]);
        if (mxGetClassID(prhs[nrhs-1]) != mxUINT32_CLASS || mxIsComplex(prhs[nrhs-1]) || (ss != 1 && ss != N))
            mexErrMsgTxt("start is of unexpected type or size");
        ss = ss == 1 ? 0 : -1;
        start = (const uint32_t *)mxGetData(prhs[nrhs-1]);
    }
	
	// Create the output array
	plhs[0] = mxCreateNumericMatrix(1, N, mxUINT32_CLASS, mxREAL);
	uint32_t *out = (uint32_t *)mxGetData(plhs[0]);
    
    // Get pointer to the input array
    const void *A = mxGetData(prhs[0]);
    
	// Call the first wrapper function according to the input image type
	switch (in_class) {
		case mxDOUBLE_CLASS:
			wrapper_func((const double *)A, op, (const double *)val, vs, start, ss, M, N, out);
			break;
		case mxSINGLE_CLASS:
			wrapper_func((const float *)A, op, (const float *)val, vs, start, ss, M, N, out);
			break;
		case mxINT8_CLASS:
			wrapper_func((const int8_t *)A, op, (const int8_t *)val, vs, start, ss, M, N, out);
			break;
		case mxUINT8_CLASS:
			wrapper_func((const uint8_t *)A, op, (const uint8_t *)val, vs, start, ss, M, N, out);
			break;
		case mxINT16_CLASS:
			wrapper_func((const int16_t *)A, op, (const int16_t *)val, vs, start, ss, M, N, out);
			break;
		case mxUINT16_CLASS:
			wrapper_func((const uint16_t *)A, op, (const uint16_t *)val, vs, start, ss, M, N, out);
			break;
		case mxINT32_CLASS:
			wrapper_func((const int32_t *)A, op, (const int32_t *)val, vs, start, ss, M, N, out);
			break;
		case mxUINT32_CLASS:
			wrapper_func((const uint32_t *)A, op, (const uint32_t *)val, vs, start, ss, M, N, out);
			break;
		case mxLOGICAL_CLASS:
			wrapper_func((const mxLogical *)A, op, (const mxLogical *)val, vs, start, ss, M, N, out);
			break;
		default:
			mexErrMsgTxt("A is of an unsupported type");
			break;
	}
	return;
}