#include "mex.h"
#include <memory.h>

// J = fill_triangles(I, P, zbuff)

// Ported to Matlab by M. Everingham

#define ALLOW_SINGLES 1

template<class T> static void FillTriangle(const T *xyv, T *img, T *tmp_buf, int ih, int iw, int id, int zbuff);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const mwSize    *dim;
	int				id, iw, ih, ntri, i;
	int 			zbuff;
    mxArray         *out_array;

	if (nrhs != 3)
		mexErrMsgTxt("3 input arguments expected.");

	if (nlhs > 1)
		mexErrMsgTxt("0 or 1 output argument expected.");

	if (ALLOW_SINGLES) {
		if ((!mxIsDouble(prhs[0]) && !mxIsSingle(prhs[0])) || mxGetClassID(prhs[0]) != mxGetClassID(prhs[1]))
			mexErrMsgTxt("inpus 1 (I) and 2 (P) must be the same type, either single or double arrays");
	} else {
		if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
			mexErrMsgTxt("inpus 1 (I) and 2 (P) must both be double arrays");
	}

	if (mxIsComplex(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) < 2 || mxGetNumberOfDimensions(prhs[0]) > 3)
		mexErrMsgTxt("input 1 (I) must be a real matrix or 3-D array");

	if (mxGetNumberOfDimensions(prhs[0]) == 2)
	{
		id = 1;
		ih = mxGetM(prhs[0]);
		iw = mxGetN(prhs[0]);
	}
	else
	{
		dim = mxGetDimensions(prhs[0]);
		id = dim[0];
		ih = dim[1];
		iw = dim[2];
	}

	if (mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1]) % (3 * (id + 2)))
		mexErrMsgTxt("input 2 (P) must be a real array with a multiple of 3*(size(I,1)+2) elements");

	ntri = mxGetNumberOfElements(prhs[1]) / (3 * (id + 2));

	if (mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1)
		mexErrMsgTxt("input 3 (zbuff) must be scalar");

	zbuff = (int)mxGetScalar(prhs[2]);

	if (zbuff != 0 && zbuff != -1 && zbuff != 1)
		mexErrMsgTxt("input 3 (zbuff) must have value -1/0/+1");
	
    if (nlhs > 0) {
        plhs[0] = mxDuplicateArray(prhs[0]);
        out_array = plhs[0];
    } else {
        // No output, so write directly into the input array. Naughty!!!!!!
        out_array = (mxArray *)prhs[0];
    }

	if (ALLOW_SINGLES && mxIsSingle(prhs[0])) {
		const float *xyv = (const float *)mxGetData(prhs[1]);
		float *out = (float *)mxGetData(out_array);
		float *tmp_buf = (float *)mxMalloc(sizeof(float)*5*id);

		for (i = 0; i < ntri; ++i)		
			FillTriangle(xyv + i * 3 * (id + 2), out, tmp_buf, iw, ih, id, zbuff);

		mxFree(tmp_buf);
	} else {
		const double *xyv = (const double *)mxGetData(prhs[1]);
		double *out = (double *)mxGetData(out_array);
		double *tmp_buf = (double *)mxMalloc(sizeof(double)*5*id);

		for (i = 0; i < ntri; ++i)		
			FillTriangle(xyv + i * 3 * (id + 2), out, tmp_buf, iw, ih, id, zbuff);

		mxFree(tmp_buf);
	}
}

template<class T> static void FillTriangle(const T *xyv, T *img, T *tmp_buf, int ih, int iw, int id, int zbuff)
{
	// note x, y etc. transposed

	const T	*v1, *v2, *v3, *vt;

	T		x_top, x_bottom, x_middle, ex;
	T		y_top, y_bottom, y_middle, ey, emy;
	T		xb_minus_xt, yb_minus_yt, ym_minus_yt, xm_minus_xt;
	T		left, left_delta, right, right_delta, width;
	T		t, t_to_b_inv_slope, t_to_m_inv_slope, m_to_b_inv_slope;
	T		recip_yb_minus_yt, recip_yb_minus_ym, recip_ym_minus_yt;
	int		iy, iy_end, iy_mid_end, ix, ix_end;
	T		*v, *v_start, *v_middle, *dvdx, *dvds;
	T		*ip, *ip_end, *ip_y;
	int		i;

	v = &tmp_buf[0*id];
	v_start = &tmp_buf[1*id];
	v_middle = &tmp_buf[2*id];
	dvdx = &tmp_buf[3*id];
	dvds = &tmp_buf[4*id];

	// sort into top, middle and bottom vertices

	v1 = &xyv[0 * (id + 2)];
	v2 = &xyv[1 * (id + 2)];
	v3 = &xyv[2 * (id + 2)];

	if (v1[0] > v2[0]) { vt = v1; v1 = v2; v2 = vt; }
	if (v1[0] > v3[0]) { vt = v1; v1 = v3; v3 = vt; }
	if (v2[0] > v3[0]) { vt = v2; v2 = v3; v3 = vt; }

	x_top = v1[1] + T(0.5);
	y_top = v1[0] + T(0.5);
	x_middle = v2[1] + T(0.5);
	y_middle = v2[0] + T(0.5);
	x_bottom = v3[1] + T(0.5);
	y_bottom = v3[0] + T(0.5);
	xb_minus_xt = x_bottom - x_top;
	yb_minus_yt = y_bottom - y_top;
	recip_yb_minus_yt = T(1.0) / yb_minus_yt;
	recip_yb_minus_ym = T(1.0) / (y_bottom - y_middle);
	ym_minus_yt = y_middle - y_top;
	recip_ym_minus_yt = T(1.0) / ym_minus_yt;
	xm_minus_xt = x_middle - x_top;

	iy = (int) (y_top + T(0.5));	// first scanline inside
	if (iy < ih)
	{
		if (iy < 1)
			iy = 1;

		ey = static_cast<T>(iy) + T(0.5) - y_top;	// vertical error from top to pixel center

		iy_end = (int) (y_bottom + T(0.5));	// first scanline outside
		if (iy_end > ih + 1)
			iy_end = ih + 1;

		iy_mid_end = (int) (y_middle + T(0.5));	// first scanline outside top "half"
		if (iy_mid_end > ih + 1)
			iy_mid_end = ih + 1;

		emy = static_cast<T>(iy_mid_end) + T(0.5) - y_middle;	// vertical error from middle to pixel center in bottom "half"

		ip_y = img + (iy - 1) * iw * id;
		t_to_b_inv_slope = recip_yb_minus_yt * xb_minus_xt;
		m_to_b_inv_slope = recip_yb_minus_ym * (x_bottom - x_middle);

		if (y_top != y_middle)
		{
			// no horizontal top

			t_to_m_inv_slope = recip_ym_minus_yt * xm_minus_xt;
			t = ym_minus_yt * recip_yb_minus_yt;
			width = xm_minus_xt - xb_minus_xt * t;

			for (i = 0; i < id; ++i)
				v_middle[i] = (v3[2 + i] - v1[2 + i]) * t + v1[2 + i];

			if (t_to_m_inv_slope < t_to_b_inv_slope)
			{
				// middle on left

				left_delta = t_to_m_inv_slope;
				right_delta = t_to_b_inv_slope;

				for (i = 0; i < id; ++i)
					dvds[i] = (v2[2 + i] - v1[2 + i]) * recip_ym_minus_yt;
			}
			else
			{
				right_delta = t_to_m_inv_slope;
				left_delta = t_to_b_inv_slope;

				for (i = 0; i < id; ++i)
					dvds[i] = (v3[2 + i] - v1[2 + i]) * recip_yb_minus_yt;
			}

			for (i = 0; i < id; ++i)
				dvdx[i] = (v2[2 + i] - v_middle[i]) / width;

			left = left_delta * ey + x_top; 
			right = right_delta * ey + x_top;

			for (i = 0; i < id; ++i)
				v_start[i] = v1[2 + i] + ey * dvds[i];

			while (iy < iy_mid_end)
			{
				ix = (int) (left + T(0.5));
				if (ix <= iw)
				{
					if (ix < 1)
						ix = 1;

					ex = static_cast<T>(ix) + T(0.5) - left;	// horizontal error from edge to pixel center

					for (i = 0; i < id; ++i)
						v[i] = v_start[i] + dvdx[i] * ex;

					ix_end = (int) (right + T(0.5));
					if (ix_end > iw + 1)
						ix_end = iw + 1;

					ip = ip_y + (ix - 1) * id;
					ip_end = ip_y + (ix_end - 1) * id;

					if (zbuff > 0)
					{
						while (ip < ip_end) 
						{
							if (v[0] < *ip)
							{
								for (i = 0; i < id; ++i)
									*(ip++) = v[i];
							}
							else
								ip += id;

							for (i = 0; i < id; ++i)
								v[i] += dvdx[i];
						}
					}
					else
						if (zbuff < 0)
						{
							while (ip < ip_end) 
							{
								if (v[0] > *ip)
								{
									for (i = 0; i < id; ++i)
										*(ip++) = v[i];
								}
								else
									ip += id;

								for (i = 0; i < id; ++i)
									v[i] += dvdx[i];
							}
						}
						else
						{
							while (ip < ip_end) 
							{
								for (i = 0; i < id; ++i)
									*(ip++) = v[i];

								for (i = 0; i < id; ++i)
									v[i] += dvdx[i];
							}
						}
				}

				for (i = 0; i < id; ++i)
					v_start[i] += dvds[i];

				left += left_delta;
				right += right_delta;
				ip_y += iw * id;
				iy++;
			} 

			// set up deltas for botthom "half" and reset x

			if (t_to_m_inv_slope < t_to_b_inv_slope)
			{
				// middle on left

				left_delta = m_to_b_inv_slope;
				left = x_middle + emy * left_delta;

				for (i = 0; i < id; ++i)
				{
					dvds[i] = (v3[2 + i] - v2[2 + i]) * recip_yb_minus_ym;
					v_start[i] = v2[2 + i] + emy * dvds[i]; 
				}
			} 
			else
			{
				right_delta = m_to_b_inv_slope; 
				right = x_middle + emy * right_delta;
			}
		}
		else 
		{
			// horizontal top

			if (x_middle < x_top)
			{
				// middle on left

				left_delta = m_to_b_inv_slope; 
				right_delta = t_to_b_inv_slope;  
				left = x_middle + emy * left_delta;
				right = x_top + emy * right_delta;  

				for (i = 0; i < id; ++i)
				{
					dvds[i] = (v3[2 + i] - v2[2 + i]) * recip_yb_minus_ym;
					v_start[i] = v2[2 + i] + emy * dvds[i];
				}
			} 
			else
			{
				right_delta = m_to_b_inv_slope; 
				left_delta = t_to_b_inv_slope;
				right = x_middle + emy * right_delta;
				left = x_top + emy * left_delta;

				for (i = 0; i < id; ++i)
				{
					dvds[i] = (v3[2 + i] - v1[2 + i]) * recip_yb_minus_yt;
					v_start[i] = v1[2 + i] + emy * dvds[i];
				}
			}
			for (i = 0; i < id; ++i)
				dvdx[i] = (v2[2 + i] - v1[2 + i]) / xm_minus_xt;
		} 

		while (iy < iy_end)
		{
			ix = (int) (left + T(0.5));
			if (ix <= iw)
			{
				if (ix < 1)
					ix = 1;

				ex = static_cast<T>(ix) + T(0.5) - left; // horizontal error from edge to pixel center

				for (i = 0; i < id; i++)
					v[i] = v_start[i] + dvdx[i] * ex;

				ix_end = (int) (right + T(0.5));
				if (ix_end > iw + 1)
					ix_end = iw + 1;

				ip = ip_y + (ix - 1) * id; 
				ip_end = ip_y + (ix_end - 1) * id;

				if (zbuff > 0)
				{
					while (ip < ip_end)
					{
						if (v[0] < *ip)
						{
							for (i = 0; i < id; ++i)
								*(ip++) = v[i];
						}
						else
							ip += id;

						for (i = 0; i < id; ++i)
							v[i] += dvdx[i];
					}
				}
				else
					if (zbuff < 0)
					{
						while (ip < ip_end)
						{
							if (v[0] > *ip)
							{
								for (i = 0; i < id; ++i)
									*(ip++) = v[i];
							}
							else
								ip += id;

							for (i = 0; i < id; ++i)
								v[i] += dvdx[i];
						}
					}
					else
					{
						while (ip < ip_end)
						{
							for (i = 0; i < id; ++i)
								*(ip++) = v[i];

							for (i = 0; i < id; ++i)
								v[i] += dvdx[i];
						}
					}
			}

			for (i = 0; i < id; ++i)
				v_start[i] += dvds[i];

			left += left_delta;
			right += right_delta;
			ip_y += iw * id;

			++iy;
		}
	}
} 
