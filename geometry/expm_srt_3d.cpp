/*
%EXPM_SRT_3D Compute a transformation matrix, given the Lie algebra vector
%
%   M = expm_srt_3d(X)
%
% Computes the transformation matrix defined by a Lie vector consisting of
% a rotation, translation and uniform scaling.
%
% The computation is done in closed form, using the formulae given in the
% paper:
% "Distances and Means of Direct Similarities"
% M-T Pham et al.
%
%IN:
%   X - DxN array, where each column specifies a different
%       transformation. X(1:3,:) are the rotation components, X(4:6,:) are
%       the translation components and X(7,:) are the scale components. D
%       can be 3 (rotation only, 6 (rotation and translation), or 7
%       (rotation, translation and scale).
%
%OUT:
%   M - 3x(3+D~=3)xN array of transformation matrices.
*/

#include <cmath>

// out - pointer to a column major matrix of size 3x4, or 3x3 if t==NULL.
// r - pointer to a 3 element vector of rotation components (x, y, z).
// t - pointer to a 3 element vector of translation components (x, y, z). Zero translation assumed if NULL.
// s - a uniform scaling component.
template<class REAL> void expm_srt_3d(REAL *out, const REAL *r, const REAL *t=NULL, const REAL s=0)
{
    // First compute the rotation part
    // Angle of rotation
    REAL theta2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    REAL theta = std::sqrt(theta2);
    REAL cos_theta = std::cos(theta);
    REAL sin_theta = std::sin(theta);
    REAL cosf = REAL(0.5);
    REAL sinc = REAL(1.0);
    if (theta2 > REAL(2.23e-16)) {
        cosf = (REAL(1.0) - cos_theta) / theta2;
        sinc = sin_theta / theta;
    }
    REAL w2[6];
    // Diagonals
    w2[0] = r[0] * r[0] - theta2;
    w2[1] = r[1] * r[1] - theta2;
    w2[2] = r[2] * r[2] - theta2;
    out[0] = w2[0] * cosf + REAL(1.0);
    out[4] = w2[1] * cosf + REAL(1.0);
    out[8] = w2[2] * cosf + REAL(1.0);
    // Off diagonals
    w2[3] = r[0] * r[1];
    out[1] = w2[3] * cosf;
    out[3] = out[1] - sinc * r[2];
    out[1] += sinc * r[2];
    w2[4] = r[0] * r[2];
    out[2] = w2[4] * cosf;
    out[6] = out[2] + sinc * r[1];
    out[2] -= sinc * r[1];
    w2[5] = r[1] * r[2];
    out[5] = w2[5] * cosf;
    out[7] = out[5] - sinc * r[0];
    out[5] += sinc * r[0];
    
    // Multiply by scale
    REAL exp_s = REAL(1.0);
    if (s) {
        exp_s = std::exp(s);
        for (int a = 0; a < 9; ++a)
            out[a] = exp_s * out[a];
    }
    
    // Now compute the translation part
    if (t) {
        REAL exp_sa = std::abs(s) < REAL(1.0e-10) ? REAL(1.0) : (exp_s - REAL(1.0)) / s;
        REAL eta_r, eta_i;
        theta2 += REAL(1e-38); // Avoid divides by 0
        if (s) {
            REAL x = exp_s * cos_theta - REAL(1.0);
            REAL y = exp_s * sin_theta;
            eta_r = exp_sa - (s * x + theta * y) / (s * s + theta2);
            eta_i = (s * y - theta * x) / (theta * (s * s + theta2) + REAL(1e-38));
        } else {
            eta_r = REAL(1.0) - sinc;
            eta_i = cosf;
        }
        eta_r /= theta2;
        
        out[9] = t[0] * (eta_r * w2[0] + exp_sa);
        out[9] += t[1] * (eta_r * w2[3] - eta_i * r[2]);
        out[9] += t[2] * (eta_r * w2[4] + eta_i * r[1]);
        out[10] = t[0] * (eta_r * w2[3] + eta_i * r[2]);
        out[10] += t[1] * (eta_r * w2[1] + exp_sa);
        out[10] += t[2] * (eta_r * w2[5] - eta_i * r[0]);
        out[11] = t[0] * (eta_r * w2[4] - eta_i * r[1]);
        out[11] += t[1] * (eta_r * w2[5] + eta_i * r[0]);
        out[11] += t[2] * (eta_r * w2[2] + exp_sa);
    }
}

#ifdef _MATLAB_
#include <mex.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#if _MATLAB_ < 805 // R2015a
extern "C" mxArray *mxCreateUninitNumericArray(mwSize ndim, const size_t *dims, mxClassID classid, mxComplexity ComplexFlag);
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs != 1)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Unexpected number of output arguments.");

	// Check argument types are valid
	if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS || mxIsComplex(prhs[0]))
		mexErrMsgTxt("Input must be a real double matrix");
    int nparams = mxGetM(prhs[0]);
    if (nparams != 3 && nparams != 6 && nparams != 7)
        mexErrMsgTxt("Input must have 3, 6 or 7 rows");
    
    // Construct the output array
    int ndims = mxGetNumberOfDimensions(prhs[0]);
	if (ndims > 19)
		mexErrMsgTxt("Input has an unsupported number of dimensions");
	size_t out_dims[20];
	out_dims[0] = 3;
	out_dims[1] = 3 + (nparams != 3);
	const int *dims = mxGetDimensions(prhs[0]);
	int nmatrices = 1;
	for (int a = 1; a < ndims; ++a) {
		out_dims[a+1] = dims[a];
		nmatrices *= dims[a];
	}
    plhs[0] = mxCreateUninitNumericArray(ndims+1, out_dims, mxDOUBLE_CLASS, mxREAL);
    double *out = (double *)mxGetData(plhs[0]);
    
    // Get the data pointers and strides
    double *r = (double *)mxGetData(prhs[0]);
    double *t = NULL;
    int t_stride = 0;
    if (nparams > 3) {
        t = r + 3;
        t_stride = nparams;
    }
    double s_ = 0;
    double *s = &s_;
    int s_stride = 0;
    if (nparams == 7) { 
        s = r + 6;
        s_stride = 7;
    }
    int out_stride = 3 * (3 + (nparams != 3));
	    
	// Call the expm function
    int a;
#pragma omp parallel for if (nmatrices > 1000) num_threads(omp_get_num_procs()) default(shared) private(a)
    for (a = 0; a < nmatrices; ++a)
        expm_srt_3d(&out[out_stride*a], &r[nparams*a], &t[t_stride*a], s[s_stride*a]);

	return;
}
#endif