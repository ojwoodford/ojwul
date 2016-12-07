#ifndef __INTERP2_METHODS_H__
#define __INTERP2_METHODS_H__

// These methods all use the following conventions:
//    1. The centre of the top left pixel is at 0,0
//    2. The image is stored in memory row-major.

// Function for correct rounding
// Add these to use numeric_limits class
#include <limits>
template<class T, class U, class V> class IM_BASE
{
public:    
	// Image dimensions
	inline int Width() { return width; }
	inline int Height() { return height; }
	inline int Channels() { return nchannels; }
    
protected:    
    const T *im;
    V *im_copy;
    int height;
    int width;
    int nchannels;
    int width_pitch;
    int plane_pitch;
    U oobv;
    V dw;
    V dh;
};

template <typename U, typename T>
static inline U saturate_cast(T val)
{
    if (std::numeric_limits<U>::is_integer && !std::numeric_limits<T>::is_integer) {
        if (std::numeric_limits<U>::is_signed)
            return val > 0 ? (val > (T)std::numeric_limits<U>::max() ? std::numeric_limits<U>::max() : static_cast<U>(val + 0.5)) : (val < (T)std::numeric_limits<U>::min() ? std::numeric_limits<U>::min() : static_cast<U>(val - 0.5));
        else
            return val > 0 ? (val > (T)std::numeric_limits<U>::max() ? std::numeric_limits<U>::max() : static_cast<U>(val + 0.5)) : 0;
    }
    return static_cast<U>(val);
}

template<class T, class U, class V> class IM_NEAR : public IM_BASE<T, U, V>
{
public:
    // Constructor
    IM_NEAR(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) {
        IM_BASE<T, U, V>::im = im_;
        IM_BASE<T, U, V>::oobv = o;
        IM_BASE<T, U, V>::width = w;
        IM_BASE<T, U, V>::dw = (V)(w) - 0.5;
        IM_BASE<T, U, V>::height = h;
        IM_BASE<T, U, V>::dh = (V)(h) - 0.5;
        IM_BASE<T, U, V>::nchannels = c;
        IM_BASE<T, U, V>::width_pitch = wp == 0 ? IM_BASE<T, U, V>::width : wp;
        IM_BASE<T, U, V>::plane_pitch = pp == 0 ? IM_BASE<T, U, V>::width_pitch * IM_BASE<T, U, V>::height : pp;
    }
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= -0.5 && X < IM_BASE<T, U, V>::dw && Y >= -0.5 && Y < IM_BASE<T, U, V>::dh) {
			// Find nearest neighbour
			int k = int(X+0.5) + int(Y+0.5) * IM_BASE<T, U, V>::width_pitch;
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, k += IM_BASE<T, U, V>::plane_pitch)
				*B = saturate_cast<U, T>(IM_BASE<T, U, V>::im[k]);
		} else {
			// Out of bounds
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch)
				*B = IM_BASE<T, U, V>::oobv;
		}
    }
    
    inline void lookup_grad(U *B, U *G, const V X, const V Y, const int out_pitch=1, const int grad_pitch=2) {}
};

template<class T, class U, class V> class IM_LIN : public IM_BASE<T, U, V>
{
public:
    // Constructor
    IM_LIN(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) {
        IM_BASE<T, U, V>::im = im_;
        IM_BASE<T, U, V>::oobv = o;
        IM_BASE<T, U, V>::width = w;
        IM_BASE<T, U, V>::dw = (V)(w) - 1.0;
        IM_BASE<T, U, V>::height = h;
        IM_BASE<T, U, V>::dh = (V)(h) - 1.0;
        IM_BASE<T, U, V>::nchannels = c;
        IM_BASE<T, U, V>::width_pitch = wp == 0 ? IM_BASE<T, U, V>::width : wp;
        IM_BASE<T, U, V>::plane_pitch = pp == 0 ? IM_BASE<T, U, V>::width_pitch * IM_BASE<T, U, V>::height : pp;
    }
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= 0.0 && Y >= 0.0 && X <= IM_BASE<T, U, V>::dw && Y <= IM_BASE<T, U, V>::dh) {
            // Compute integer coordinates and offsets
            int y = (int)(Y);
            int x = (int)(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == IM_BASE<T, U, V>::dw) {
                x = x > 0 ? --x : 0;
                u = V(1.0);
            }
            if (Y == IM_BASE<T, U, V>::dh) {
                y = y > 0 ? --y : 0;
                v = V(1.0);
            }
            // Compute the linear index
            int k = x + y * IM_BASE<T, U, V>::width_pitch;
            // For each image channel...
            for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, k += IM_BASE<T, U, V>::plane_pitch) {
                // Do the interpolation
                V out = IM_BASE<T, U, V>::im[k] + (IM_BASE<T, U, V>::im[k+1] - IM_BASE<T, U, V>::im[k]) * u;
                out += ((IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch] - out) + (IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch+1] - IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch]) * u) * v;
                *B = saturate_cast<U, V>(out);
            }
		} else {
			// Out of bounds
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch)
				*B = IM_BASE<T, U, V>::oobv;
		}
    }
    
    // Lookup value and gradient function
    inline void lookup_grad(U *B, U *G, const V X, const V Y, const int out_pitch=1, const int grad_pitch=2) {
        if (X >= 0.0 && Y >= 0.0 && X <= IM_BASE<T, U, V>::dw && Y <= IM_BASE<T, U, V>::dh) {
            // Compute integer coordinates and offsets
            int y = (int)(Y);
            int x = (int)(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == IM_BASE<T, U, V>::dw) {
                x = x > 0 ? --x : 0;
                u = V(1.0);
            }
            if (Y == IM_BASE<T, U, V>::dh) {
                y = y > 0 ? --y : 0;
                v = V(1.0);
            }
            // Compute the linear index
            int k = x + y * IM_BASE<T, U, V>::width_pitch;
            // For each image channel...
            for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, G += grad_pitch, k += IM_BASE<T, U, V>::plane_pitch) {
                // Sample the image
                V v00 = (V)IM_BASE<T, U, V>::im[k];
                V v10 = (V)IM_BASE<T, U, V>::im[k+1];
                V v01 = (V)IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch];
                V v11 = (V)IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch+1];
                // Linearly interpolate
                V d0 = v10 - v00;
                V d1 = v11 - v01;
                V out = v00 + d0 * u;
                out += ((v01 - out) + d1 * u) * v;
                *B = saturate_cast<U, V>(out);
                out = d0 + (d1 - d0) * v;
                G[1] = saturate_cast<U, V>(out);
                d0 = v01 - v00;
                d1 = v11 - v10;
                out = d0 + (d1 - d0) * u;
                G[0] = saturate_cast<U, V>(out);
            }
		} else {
			// Out of bounds
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, G += grad_pitch) {
                *B = IM_BASE<T, U, V>::oobv;
                G[0] = U(0);
                G[1] = U(0);
            }
		}
    }
};

// This shifted interpolation method is described in:
// "Linear Interpolation Revitalized", T. Blu, P. Thevenaz, M. Unser.
// IEEE Transactions on Image Processing, vol. 13, no. 5, May 2004.
#define __SHIFT_FACTOR__ 0.21132486540518713 // 0.5 * (1 - sqrt(3) / 3)
template<class T, class U, class V> class IM_SHIFT : public IM_BASE<T, U, V>
{
public:
    // Constructor
    IM_SHIFT(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) {
        IM_BASE<T, U, V>::oobv = o;
        IM_BASE<T, U, V>::width = w;
        IM_BASE<T, U, V>::dw = (V)(w) - 1.0;
        IM_BASE<T, U, V>::height = h;
        IM_BASE<T, U, V>::dh = (V)(h) - 1.0;
        IM_BASE<T, U, V>::nchannels = c;
        wp = wp == 0 ? IM_BASE<T, U, V>::width : wp;
        pp = pp == 0 ? wp * IM_BASE<T, U, V>::height : pp;
        IM_BASE<T, U, V>::width_pitch = IM_BASE<T, U, V>::width;
        IM_BASE<T, U, V>::plane_pitch = IM_BASE<T, U, V>::width_pitch * IM_BASE<T, U, V>::height;
        
        // Allocate memory for the new coefficients
        IM_BASE<T, U, V>::im_copy = new V[IM_BASE<T, U, V>::plane_pitch*IM_BASE<T, U, V>::nchannels];
        
        const V factor = 1.0 / (1.0 - __SHIFT_FACTOR__);
        const V pole = __SHIFT_FACTOR__ / (__SHIFT_FACTOR__ - 1.0);
        
        // Conversion along rows
        for (c = 0; c < IM_BASE<T, U, V>::nchannels; ++c) {
            for (h = 0; h < IM_BASE<T, U, V>::height; ++h) {
                V *im_ptr = &IM_BASE<T, U, V>::im_copy[c*IM_BASE<T, U, V>::plane_pitch+h*IM_BASE<T, U, V>::width_pitch];
                const T *im__ptr = &im_[c*pp+h*wp];
                im_ptr[0] = (V)(im__ptr[0]);
                for (w = 1; w < IM_BASE<T, U, V>::width; ++w)
                    im_ptr[w] = pole * im_ptr[w-1] + factor * (V)(im__ptr[w]);
            }
        }
        
        // In place conversion along columns
        for (c = 0; c < IM_BASE<T, U, V>::nchannels; ++c) {
            for (w = 0; w < IM_BASE<T, U, V>::width; ++w) {
                V *im_ptr = &IM_BASE<T, U, V>::im_copy[c*IM_BASE<T, U, V>::plane_pitch+w];
                for (h = 1; h < IM_BASE<T, U, V>::height; ++h, im_ptr += IM_BASE<T, U, V>::width_pitch)
                    im_ptr[IM_BASE<T, U, V>::width_pitch] = pole * im_ptr[0] + factor * im_ptr[IM_BASE<T, U, V>::width_pitch];
            }
        }
    }
    // Destructor
    ~IM_SHIFT() {
        delete [] IM_BASE<T, U, V>::im_copy;
    }
    
    // Lookup function
    inline void lookup(U *B, V X, V Y, const int out_pitch=1) {
        // Do the shift
        X -= __SHIFT_FACTOR__;
        Y -= __SHIFT_FACTOR__;
        // Check in bounds
        if (X >= 0.0 && Y >= 0.0 && X <= IM_BASE<T, U, V>::dw && Y <= IM_BASE<T, U, V>::dh) {
            // Compute integer coordinates and offsets
            int y = (int)(Y);
            int x = (int)(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == IM_BASE<T, U, V>::dw) {
                x = x > 0 ? --x : 0;
                u = V(1.0);
            }
            if (Y == IM_BASE<T, U, V>::dh) {
                y = y > 0 ? --y : 0;
                v = V(1.0);
            }
            // Compute the linear index
            int k = x + y * IM_BASE<T, U, V>::width_pitch;
            // For each image channel...
            for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, k += IM_BASE<T, U, V>::plane_pitch) {
                // Do the interpolation
                V out = IM_BASE<T, U, V>::im_copy[k] + (IM_BASE<T, U, V>::im_copy[k+1] - IM_BASE<T, U, V>::im_copy[k]) * u;
                out += ((IM_BASE<T, U, V>::im_copy[k+IM_BASE<T, U, V>::width_pitch] - out) + (IM_BASE<T, U, V>::im_copy[k+IM_BASE<T, U, V>::width_pitch+1] - IM_BASE<T, U, V>::im_copy[k+IM_BASE<T, U, V>::width_pitch]) * u) * v;
                *B = saturate_cast<U, V>(out);
            }
		} else {
			// Out of bounds
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch)
				*B = IM_BASE<T, U, V>::oobv;
		}
    }
    
    // Lookup value and gradient function
    inline void lookup_grad(U *B, U *G, V X, V Y, const int out_pitch=1, const int grad_pitch=2) {
        // Do the shift
        X -= __SHIFT_FACTOR__;
        Y -= __SHIFT_FACTOR__;
        if (X >= 0.0 && Y >= 0.0 && X <= IM_BASE<T, U, V>::dw && Y <= IM_BASE<T, U, V>::dh) {
            // Compute integer coordinates and offsets
            int y = (int)(Y);
            int x = (int)(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == IM_BASE<T, U, V>::dw) {
                x = x > 0 ? --x : 0;
                u = V(1.0);
            }
            if (Y == IM_BASE<T, U, V>::dh) {
                y = y > 0 ? --y : 0;
                v = V(1.0);
            }
            // Compute the linear index
            int k = x + y * IM_BASE<T, U, V>::width_pitch;
            // For each image channel...
            for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, G += grad_pitch, k += IM_BASE<T, U, V>::plane_pitch) {
                // Sample the image
                V v00 = IM_BASE<T, U, V>::im_copy[k];
                V v10 = IM_BASE<T, U, V>::im_copy[k+1];
                V v01 = IM_BASE<T, U, V>::im_copy[k+IM_BASE<T, U, V>::width_pitch];
                V v11 = IM_BASE<T, U, V>::im_copy[k+IM_BASE<T, U, V>::width_pitch+1];
                // Linearly interpolate
                V d0 = v10 - v00;
                V d1 = v11 - v01;
                V out = v00 + d0 * u;
                out += ((v01 - out) + d1 * u) * v;
                *B = saturate_cast<U, V>(out);
                G[1] = saturate_cast<U, V>(d0 + (d1 - d0) * v);
                d0 = v01 - v00;
                d1 = v11 - v10;
                G[0] = saturate_cast<U, V>(d0 + (d1 - d0) * u);
            }
		} else {
			// Out of bounds
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, G += grad_pitch) {
                *B = IM_BASE<T, U, V>::oobv;
                G[0] = U(0);
                G[1] = U(0);
            }
		}
    }
};

template<class T, class U, class V> class IM_CUB : public IM_BASE<T, U, V>
{
public:
    // Constructor
    IM_CUB(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) {
        IM_BASE<T, U, V>::im = im_;
        IM_BASE<T, U, V>::oobv = o;
        IM_BASE<T, U, V>::width = w;
        IM_BASE<T, U, V>::dw = (V)(w) - 2.0;
        IM_BASE<T, U, V>::height = h;
        IM_BASE<T, U, V>::dh = (V)(h) - 2.0;
        IM_BASE<T, U, V>::nchannels = c;
        IM_BASE<T, U, V>::width_pitch = wp == 0 ? IM_BASE<T, U, V>::width : wp;
        IM_BASE<T, U, V>::plane_pitch = pp == 0 ? IM_BASE<T, U, V>::width_pitch * IM_BASE<T, U, V>::height : pp;
    }
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= 2.0 && X < IM_BASE<T, U, V>::dw && Y >= 2.0 && Y < IM_BASE<T, U, V>::dh) {
            // Bicubicly interpolate
            V b[4], d[4], u[3], v[3];
			int x = (int)X;
			int y = (int)Y;
			u[0] = X - x;
			v[0] = Y - y;
			u[1] = u[0] * u[0];
			v[1] = v[0] * v[0];
			u[2] = u[1] * u[0];
			v[2] = v[1] * v[0];
			int k = x - 1 + (y - 1) * IM_BASE<T, U, V>::width_pitch;
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, k += IM_BASE<T, U, V>::plane_pitch) {
				V a;
                for (int m = 0, n = k; m < 4; ++m, n += IM_BASE<T, U, V>::width_pitch) {
					d[0] = (V)IM_BASE<T, U, V>::im[n+0];
					d[1] = (V)IM_BASE<T, U, V>::im[n+1];
					d[2] = (V)IM_BASE<T, U, V>::im[n+2];
					d[3] = (V)IM_BASE<T, U, V>::im[n+3];
					a = (d[3] + d[1]) - (d[2] + d[0]);
					b[m] = u[2] * a + u[1] * ((d[0] - d[1]) - a) + u[0] * (d[2] - d[0]) + d[1];
				}
				a = (b[3] + b[1]) - (b[2] + b[0]);
				*B = saturate_cast<U, V>(v[2] * a + v[1] * ((b[0] - b[1]) - a) + v[0] * (b[2] - b[0]) + b[1]);
			}
		} else {
			// Out of bounds
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch)
				*B = IM_BASE<T, U, V>::oobv;
		}
    }
    
    inline void lookup_grad(U *B, U *G, const V X, const V Y, const int out_pitch=1, const int grad_pitch=2) {}
};

#endif
