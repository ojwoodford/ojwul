/* [CX sse I] = fast_kmeans(X, params, [CX])

IN:
	X - MxN matrix of N input vectors of dimension M.
    params - [nclusters max_iters min_sse_delta]

OUT:
	CX - MxP matrix of updated cluster centres of X.
	sse - Sum of squared errors of distances to cluster centres.
	I - Nx1 uint32 vector of cluster indices each input vector was in this iteration

*/

#include <mex.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// Define types
#ifdef _MSC_VER
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

template<class T> static void initialize_centres(const T *X, T *CX, int ndims, int nvecs, int nclusters);
template<class T> static double wrapper_func(const T *X, T *CX, int32_t *I, int ndims, int nvecs, int nclusters, int max_iters, double min_delta, double robust_thresh);
template<class T> static inline double kmeans_func(const T *X_, T *CX, T *CX_temp, int32_t *I, int *CN, T *S, T *md, T robust_thresh, int ndims, int nvecs, int nclusters);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2 || nrhs > 3)
        mexErrMsgTxt("Unexpected number of input arguments.");
    if (nlhs < 1 || nlhs > 3)
        mexErrMsgTxt("Unsupported number output arguments expected.");
    
    if (mxIsComplex(prhs[0]) ||	mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("X must be a real matrix");
    
    mxClassID class_X = mxGetClassID(prhs[0]);
    int ndims = mxGetM(prhs[0]);
    if (ndims < 1)
        mexErrMsgTxt("size(X, 1) must be greater than 0.");
    int nvecs = mxGetN(prhs[0]);
    const void *X = mxGetData(prhs[0]);
    
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||	mxGetNumberOfElements(prhs[1]) != 4)
        mexErrMsgTxt("options must be a real 1x4 double vector");
    const double *options = mxGetPr(prhs[1]);
    int nclusters = (int)options[0];
    int max_iters = (int)options[1];
    double min_delta = options[2];
    double robust_thresh = options[3];

    // Get or create our cluster centre array
    void *CX;
    if (nrhs < 3) {
        plhs[0] = mxCreateNumericMatrix(ndims, nclusters, class_X, mxREAL);
        CX = mxGetData(plhs[0]);
        switch (class_X) {
            case mxDOUBLE_CLASS:
                initialize_centres((const double *)X, (double *)CX, ndims, nvecs, nclusters);
                break;
            case mxSINGLE_CLASS:
                initialize_centres((const float *)X, (float *)CX, ndims, nvecs, nclusters);
                break;
            case mxINT32_CLASS:
                initialize_centres((const int32_t *)X, (int32_t *)CX, ndims, nvecs, nclusters);
                break;
            default:
                mexErrMsgTxt("X is of an unsupported type");
                break;
        }
    } else {
        if (ndims != mxGetM(prhs[2]) || nclusters != mxGetN(prhs[2]) || class_X != mxGetClassID(prhs[2]) || mxIsComplex(prhs[2]))
            mexErrMsgTxt("CX must be an (ndims)x(nclusters) real matrix of the same class as X.");
        plhs[0] = mxDuplicateArray(prhs[2]);
        CX = mxGetData(plhs[0]);
    }
    
    // Create our index array
    mxArray *I_array = mxCreateNumericMatrix(1, nvecs, mxINT32_CLASS, mxREAL);
    int32_t *I = (int32_t *)mxGetData(I_array);

	double sse;
    switch (class_X) {
        case mxDOUBLE_CLASS:
            sse = wrapper_func((const double *)X, (double *)CX, I, ndims, nvecs, nclusters, max_iters, min_delta, robust_thresh);
            break;
        case mxSINGLE_CLASS:
            sse = wrapper_func((const float *)X, (float *)CX, I, ndims, nvecs, nclusters, max_iters, min_delta, robust_thresh);
            break;
        case mxINT32_CLASS:
            sse = wrapper_func((const int32_t *)X, (int32_t *)CX, I, ndims, nvecs, nclusters, max_iters, min_delta, robust_thresh);
            break;
        default:
            mexErrMsgTxt("X is of an unsupported type");
            break;
    }

    if (nlhs > 2) {
        // Update the cluster centres for Matlab notation
        for (int i = 0; i < nvecs; ++i)
            ++I[i];
        plhs[2] = I_array;
    } else {
        mxDestroyArray(I_array);
    }
    
	// Output the sse
	if (nlhs > 1)
		plhs[1] = mxCreateDoubleScalar(sse);
    
    return;
}

template<class T> static void initialize_centres(const T *X, T *CX, int ndims, int nvecs, int nclusters)
{
    // Initialise our cluster centres
    mxArray *mx_arrays[2];
    mx_arrays[1] = mxCreateDoubleScalar((double)nvecs);
    mexCallMATLAB(1, mx_arrays, 1, mx_arrays+1, "randperm");
    double *order = mxGetPr(mx_arrays[0]);
    T *V = CX;
    for (int i = 0; i < nclusters; ++i) {
        const T *V2 = &X[(((int)order[i])-1)*ndims];
        for (int j = 0; j < ndims; ++j)
            V[j] = V2[j];
        V += ndims;
    }
    mxDestroyArray(mx_arrays[0]);
    mxDestroyArray(mx_arrays[1]);
    return;
}

template<class T> static double wrapper_func(const T *X, T *CX, int32_t *I, int ndims, int nvecs, int nclusters, int max_iters, double min_delta, double robust_thresh)
{
    // Malloc our memory
    T *CX_temp = (T *)mxCalloc(ndims*nclusters, sizeof(T));
    int *CN = (int *)mxCalloc(nclusters, sizeof(int));
    T *S = (T *)mxMalloc(nclusters*sizeof(T));
    T *md = (T *)mxMalloc(nvecs*sizeof(T));
    
    // Main loop
    mxArray *mx_arrays[2];
    mx_arrays[1] = mxCreateDoubleScalar(0.0);
    mx_arrays[0] = mxCreateString("Computing k-means clusters");
    double *progress = mxGetPr(mx_arrays[1]);
    double sse_new, sse = 1.0e300;
    double recip_iters = 1.0 / (double)max_iters;
    T rt = (T)robust_thresh;
    for (int iter = 0; iter < max_iters; ++iter) {
        // Display progress
        *progress = recip_iters * (double)iter;
        mexCallMATLAB(0, NULL, 2, mx_arrays, "ojw_progressbar");
        
        // Update centres
        try {
            sse_new = kmeans_func(X, CX, CX_temp, I, CN, S, md, rt, ndims, nvecs, nclusters);
        }
        catch (char *str) {
            mexErrMsgTxt(str);
        }
        
        // Check convergence criterion
        if ((sse - sse_new) <= min_delta)
            break;
        sse = sse_new;
    }
    
    // Close the progress bar
    *progress = 1;
    mexCallMATLAB(0, NULL, 2, mx_arrays, "ojw_progressbar");
    
    // Free memory
    mxDestroyArray(mx_arrays[0]);
    mxDestroyArray(mx_arrays[1]);
    mxFree(md);
    mxFree(S);
    mxFree(CN);
    mxFree(CX_temp);
    
    return sse_new;
}
    
template<class T> static inline double kmeans_func(const T *X_, T *CX, T *CX_temp, int32_t *I, int *CN, T *S, T *md, T robust_thresh, int ndims, int nvecs, int nclusters)
{    
	int i, j, k;

    // For each cluster, compute the distance to the closest cluster
    for (i = 0; i < nclusters; ++i)
        S[i] = (T)1.0e300;
    const T *V = CX;
    for (i = 0; i < nclusters-1; ++i) {
        const T *V2 = V + ndims;
        for (j = i+1; j < nclusters; ++j) {
            T dist = 0.0;
            for (k = 0; k < ndims; ++k)
                dist += (V[k] - V2[k]) * (V[k] - V2[k]);
            if (dist < S[i])
                S[i] = dist;
            if (dist < S[j])
                S[j] = dist;
            V2 += ndims;
        }
        V += ndims;
    }
    for (i = 0; i < nclusters; ++i)
        S[i] /= 4;
    
	int ndims4 = 4 * ((ndims - 1) / 4) + 1;
    double sse = 0.0;
#pragma  omp parallel for if (nvecs > 300) num_threads(omp_get_num_procs()) default(shared) private(i,j,k,V) reduction(+:sse)
	for (k = 0; k < nvecs; ++k) {
        // Start with cluster centre from previous iteration
        T min_dist = 0.0;
        int C_min = I[k];
		const T *X = &X_[k*ndims];
        V = &CX[C_min*ndims];
        for (i = 0; i < ndims; ++i)
            min_dist += (V[i] - X[i]) * (V[i] - X[i]);
        
        // Check whether this data point could be closer to any other cluster centre
        if (min_dist > S[C_min]) {
            // Need to go through other clusters
            int Cmin_old = C_min;
            V = CX;
            for (j = 0; j < Cmin_old; ++j, V += ndims) {
                // Compute first value and check - if data is PCA'd this will speed things up
                T dist = (V[0] - X[0]) * (V[0] - X[0]);
                if (dist >= min_dist)
                    continue;
                
                // Go through the rest in 4s
                for (i = 1; i < ndims4; i += 4) {
                    dist += (V[i] - X[i]) * (V[i] - X[i]);
                    dist += (V[i+1] - X[i+1]) * (V[i+1] - X[i+1]);
                    dist += (V[i+2] - X[i+2]) * (V[i+2] - X[i+2]);
                    dist += (V[i+3] - X[i+3]) * (V[i+3] - X[i+3]);
                    if (dist >= min_dist)
                        goto skip_point1;
                }
                
                // Finish up the rest not a multiple of 4
                for (i = ndims4; i < ndims; ++i)
                    dist += (V[i] - X[i]) * (V[i] - X[i]);
                if (dist < min_dist) {
                    min_dist = dist;
                    C_min = j;
                }
skip_point1:
                continue;
            }
            
            // Leave a gap where the starting cluster is
            V += ndims;
            
            for (j = Cmin_old + 1; j < nclusters; ++j, V += ndims) {
                // Compute first value and check - if data is PCA'd this will speed things up
                T dist = (V[0] - X[0]) * (V[0] - X[0]);
                if (dist >= min_dist)
                    continue;
                
                // Go through the rest in 4s
                for (i = 1; i < ndims4; i += 4) {
                    dist += (V[i] - X[i]) * (V[i] - X[i]);
                    dist += (V[i+1] - X[i+1]) * (V[i+1] - X[i+1]);
                    dist += (V[i+2] - X[i+2]) * (V[i+2] - X[i+2]);
                    dist += (V[i+3] - X[i+3]) * (V[i+3] - X[i+3]);
                    if (dist >= min_dist)
                        goto skip_point2;
                }
                
                // Finish up the rest not a multiple of 4
                for (i = ndims4; i < ndims; ++i)
                    dist += (V[i] - X[i]) * (V[i] - X[i]);
                if (dist < min_dist) {
                    min_dist = dist;
                    C_min = j;
                }
skip_point2:
                continue;
            }
            I[k] = (int32_t)C_min;
        }
        md[k] = min_dist;
        sse += (double)min_dist;
    }

    if (robust_thresh) {
        for (k = 0; k < nvecs; ++k) {
            // Check this is within the threshold
            if (md[k] > robust_thresh)
                continue;
            // Add the vectors to the relevant cluster centre
            ++CN[I[k]];
            V = &X_[k*ndims];
            T *U = &CX_temp[I[k]*ndims];
            for (int i = 0; i < ndims; ++i)
                U[i] += V[i];
        }
    } else {
        for (k = 0; k < nvecs; ++k) {
            // Add the vectors to the relevant cluster centre
            ++CN[I[k]];
            V = &X_[k*ndims];
            T *U = &CX_temp[I[k]*ndims];
            for (int i = 0; i < ndims; ++i)
                U[i] += V[i];
        }
    }
    
	// Compute the new cluster centres
    for (int j = 0; j < nclusters; ++j, CX_temp += ndims, CX += ndims) {
        if (CN[j] > 1) {
            if ((T)0.3) {
                // Floating point aritmetic
                T recip = 1.0 / (T)CN[j];
                for (int i = 0; i < ndims; ++i) {
                    CX[i] = CX_temp[i] * recip;
                    CX_temp[i] = 0; // Clear for next time
                }
            } else {
                // Integer arithmetic
                for (int i = 0; i < ndims; ++i) {
                    CX[i] = CX_temp[i] / CN[j];
                    CX_temp[i] = 0; // Clear for next time
                }
            }
        } else if (CN[j] == 0) {
            // Empty cluster - assign to a random vector
            int pick = (int)(((double)rand()) * ((double)nvecs) / (((double)RAND_MAX) + 1.0));
            V = &X_[pick*ndims];
            for (int i = 0; i < ndims; i++)
                CX[i] = V[i];
        } else { // CN[j] == 1
            for (int i = 0; i < ndims; i++) {
                CX[i] = CX_temp[i];
                CX_temp[i] = 0; // Clear for next time
            }
        }
        // Clear for next time
        CN[j] = 0;
    }    
    return sse;
}
