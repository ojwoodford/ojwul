/* [en L] = dp(U, E); */
#include "mex.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_procs()  1
#define omp_get_thread_num() 0
#endif

/* Define types */
#include <stdint.h>

#if _MATLAB_ < 805 /* R2015a */
extern "C" mxArray *mxCreateUninitNumericMatrix(mwSize m, mwSize n, mxClassID classid, mxComplexity ComplexFlag);
#endif
 
template<class T> static T main_func(const T *U, const T *E, int32_t *L, const int nLabels, const int nNodes, const int Eoffset, int8_t *tmp_buf)
{
    T *runningE[2];
    T minE, currE;
    int *bookKeeping, *bookPtr;
    int l, l2, n, curr, minLab;
    
    /* Create temporary arrays */
    runningE[0] = (T *)tmp_buf;
    runningE[1] = runningE[0] + nLabels;
    bookKeeping = (int *)(runningE[1] + nLabels);
    
    /* Initialise the first column of energies */
    for (l = 0; l < nLabels; ++l)
        runningE[0][l] = *U++;
    
    /* Go through each following column */
    curr = 0;
    bookPtr = bookKeeping;
    for (n = 1; n < nNodes; ++n) {
        /* For each label, compute the min energy and source label */
        for (l = 0; l < nLabels; ++l) {
            minLab = 0;
            minE = runningE[curr][0] + *E++;
            for (l2 = 1; l2 < nLabels; ++l2) {
                currE = runningE[curr][l2] + *E++;
                if (currE < minE) {
                    minE = currE;
                    minLab = l2;
                }
            }
            runningE[1-curr][l] = minE + *U++;
            *bookPtr++ = minLab;
        }
        curr = 1 - curr;
        E -= Eoffset;
    }
    
    /* Find the minimum energy */
    minLab = 0;
    minE = runningE[curr][0];
    for (l = 1; l < nLabels; ++l) {
        if (runningE[curr][l] < minE) {
            minE = runningE[curr][l];
            minLab = l;
        }
    }
        
    /* Compute the optimal labelling */
    L = L + (nNodes - 1);
    *L-- = minLab + 1;
    bookPtr = bookKeeping + nLabels*(nNodes-2);
    while (bookPtr >= bookKeeping) {
        minLab = bookPtr[minLab];
        *L-- = minLab + 1;
        bookPtr -= nLabels;
    }
    
    /* Save the energy */
    return minE;
}

template<class T> static void wrapper_func(const T *U, const T *E, int32_t *L, const int nLabels, const int nNodes, const int Eoffset, T *en, const int nChains, const int Eoff) 
{
    int n, m = 2 * nLabels * sizeof(T) + nLabels * (nNodes-1) * sizeof(int);
    int8_t *tmp_buf = (int8_t *)mxMalloc(m*omp_get_num_procs());
    /* Use all the processors we can */
    if (en) {
#pragma omp parallel for if (nChains > 10) num_threads(omp_get_num_procs())
        for (n = 0; n < nChains; ++n) {
            en[n] = main_func(&U[n*nNodes*nLabels], &E[n*Eoff], &L[n*nNodes], nLabels, nNodes, Eoffset, tmp_buf+m*omp_get_thread_num());
        }
    } else {
#pragma omp parallel for if (nChains > 10) num_threads(omp_get_num_procs())
        for (n = 0; n < nChains; ++n) {
            main_func(&U[n*nNodes*nLabels], &E[n*Eoff], &L[n*nNodes], nLabels, nNodes, Eoffset, tmp_buf+m*omp_get_thread_num());
        }
    }
    mxFree(tmp_buf);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int nLabels, nNodes, nChains, a, n, Eoffset, Eoff;
    const void *U, *E;
    void *en = NULL;
    const mwSize *dims;
    int32_t *L;
    mxClassID classid;
    
    /* Check number of arguments */
    if (nrhs != 2)
        mexErrMsgTxt("Unexpected number of input arguments.");
    if (nlhs < 1 || nlhs > 2)
        mexErrMsgTxt("Unexpected number of output arguments.");
    
    /* Check input types are valid */
    classid = mxGetClassID(prhs[0]);
    if (mxGetClassID(prhs[1]) != classid)
        mexErrMsgTxt("E must be the same type as U.");
    n = mxGetNumberOfDimensions(prhs[0]);
    if (n > 3)
        mexErrMsgTxt("U has unexpected dimensions");
    dims = mxGetDimensions(prhs[0]);
    nLabels = dims[0];
    nNodes = dims[1];
    if (nNodes*nLabels < 1)
        mexErrMsgTxt("U cannot be empty.");
    if (n == 3) {
        nChains = dims[2];
    } else {
        nChains = 1;
    }
    n = mxGetNumberOfDimensions(prhs[1]);
    if (n > 4)
        mexErrMsgTxt("E has unexpected dimensions");
    dims = mxGetDimensions(prhs[1]);
    if (dims[0] != nLabels || dims[1] != nLabels)
        mexErrMsgTxt("E has unexpected dimensions");
    if (n > 2) {
        if (dims[2] != nNodes-1)
            mexErrMsgTxt("E has unexpected dimensions");
        if (n == 4) {
            if (dims[3] != nChains)
                mexErrMsgTxt("E has unexpected dimensions");
        } else {
            if (1 != nChains)
                mexErrMsgTxt("E has unexpected dimensions");
        }
        Eoffset = 0;
        Eoff = nLabels*nLabels*(nNodes-1);
    } else {
        Eoffset = nLabels*nLabels;
        Eoff = 0;
    }
    
    /* Get array pointers */
    U = (const void *)mxGetData(prhs[0]);
    E = (const void *)mxGetData(prhs[1]);
    
    /* Create the output matrices */
    plhs[0] = mxCreateUninitNumericMatrix(nNodes, nChains, mxINT32_CLASS, mxREAL);
    L = (int32_t *)mxGetData(plhs[0]);
    if (nlhs > 1) {
        plhs[1] = mxCreateUninitNumericMatrix(1, nChains, classid, mxREAL);
        en = (void *)mxGetData(plhs[1]);
    }
    
    /* Call the wrapper function according to type */
    switch (classid) {
        case mxDOUBLE_CLASS:
			wrapper_func((const double *)U, (const double *)E, L, nLabels, nNodes, Eoffset, (double *)en, nChains, Eoff);
			break;
		case mxSINGLE_CLASS:
            wrapper_func((const float *)U, (const float *)E, L, nLabels, nNodes, Eoffset, (float *)en, nChains, Eoff);
            break;
		case mxUINT32_CLASS:
            wrapper_func((const uint32_t *)U, (const uint32_t *)E, L, nLabels, nNodes, Eoffset, (uint32_t *)en, nChains, Eoff);
            break;
        default:
			mexErrMsgTxt("U and E are of an unsupported type");
			break;            
    }
}
