#ifndef __INTERP2_METHODS_H__
#define __INTERP2_METHODS_H__
#include <limits> // For std::numeric_limits
#include <utility> // For std::pair
#include <cstdlib> // For std::abs
#include <algorithm> // For std::min and std::max
#include <vector>

// These methods all use the following conventions:
//    1. The centre of the top left pixel is at 0,0
//    2. The image is stored in memory row-major.

template<class T, class U, class V> class IM_BASE
{
public:
	// Constructor
	IM_BASE(const T* im_, U o, int w, int h, int c, int wp, int pp) :
		im(im_),
		oobv(o),
		height(h),
		width(w),
		nchannels(c),
		width_pitch(wp == 0 ? w : wp),
		plane_pitch(pp == 0 ? width_pitch * h : pp) {}

    // Image dimensions
	inline int Width() const { return width; }
	inline int Height() const { return height; }
	inline int Channels() const { return nchannels; }
	inline int WidthPitch() const { return width_pitch; }
	inline int PlanePitch() const { return plane_pitch; }
    
    // Other helpers
	inline U OOBV() const { return oobv; }
	inline void SetOOBV(U new_oobv) { oobv = new_oobv; }

	inline int ind_symmetric(int x, int y) const // Lookup with symmetric padding
	{
		if (x < 0 || x >= width) {
			x = x % (2 * width);
			x = (x + 2 * width) % (2 * width);
			if (x >= width)
				x = 2 * width - 1 - x;
		}
		if (y < 0 || y >= height) {
			y = y % (2 * height);
			y = (y + 2 * height) % (2 * height);
			if (y >= height)
				y = 2 * height - 1 - y;
		}
		return x + y * width_pitch;
	}

	// Get pixel values
	inline const T& operator[](int ind) const { return im[ind]; }
	inline V at(int ind) const { return static_cast<V>(im[ind]); }
    
private:    
    const T* const im;
    const int height;
    const int width;
    const int nchannels;
    const int width_pitch;
    const int plane_pitch;
    U oobv;
};

// Function for correct rounding
template <typename U, typename T>
static inline U saturate_cast(T val)
{
    if (std::numeric_limits<U>::is_integer && !std::numeric_limits<T>::is_integer) {
        val += (val > static_cast<T>(0)) ? static_cast<T>(0.5) : -static_cast<T>(0.5);
        val = std::min(val, static_cast<T>(std::numeric_limits<U>::max()));
        val = std::max(val, static_cast<T>(std::numeric_limits<U>::min()));
    }
    return static_cast<U>(val);
}

template<class T, class U, class V> class IM_NEAR
{
public:
    // Constructor
    IM_NEAR(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) : 
		im(im_, o, w, h, c, wp, pp),
		dw(static_cast<V>(w) - static_cast<V>(0.5)),
		dh(static_cast<V>(h) - static_cast<V>(0.5)) {}
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= -0.5 && X < dw && Y >= -0.5 && Y < dh) {
			// Find nearest neighbour
			int k = static_cast<int>(X + static_cast<V>(0.5)) + static_cast<int>(Y + static_cast<V>(0.5)) * im.WidthPitch();
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch, k += im.PlanePitch())
				*B = saturate_cast<U, T>(im[k]);
		} else {
			// Out of bounds
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch)
				*B = im.OOBV();
		}
    }
    
    // Lookup value and gradient function
    inline void lookup_grad(U *B, U *G, const V X, const V Y, const int out_pitch=1, const int grad_pitch=2) {
        if (X >= -0.5 && X < dw && Y >= -0.5 && Y < dh) {
			// Find nearest neighbour
            int y = static_cast<int>(Y + static_cast<V>(0.5));
            int x = static_cast<int>(X + static_cast<V>(0.5));
			int k = x + y * im.WidthPitch();
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch, G += grad_pitch, k += im.PlanePitch())
            {
				*B = saturate_cast<U, T>(im[k]);
                if (y == 0)
                    G[0] = saturate_cast<U, V>(im.at(k+im.WidthPitch()) - im.at(k));
                else if (y == im.Height()-1)
                    G[0] = saturate_cast<U, V>(im.at(k) - im.at(k-im.WidthPitch()));
                else
                    G[0] = saturate_cast<U, V>(static_cast<V>(0.5) * (im.at(k+im.WidthPitch()) - im.at(k-im.WidthPitch())));
                if (x == 0)
                    G[1] = saturate_cast<U, V>(im.at(k+1) - im.at(k));
                else if (x == im.Width()-1)
                    G[1] = saturate_cast<U, V>(im.at(k) - im.at(k-1));
                else
                    G[1] = saturate_cast<U, V>(static_cast<V>(0.5) * (im.at(k+1) - im.at(k-1)));
            }
		} else {
			// Out of bounds
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch, G += grad_pitch) {
                *B = im.OOBV();
                G[0] = static_cast<U>(0);
                G[1] = static_cast<U>(0);
            }
		}
    }

private:
	const IM_BASE<T, U, V> im;
	const V dw;
	const V dh;
};

template<class T, class U, class V> class IM_LIN
{
public:
    // Constructor
    IM_LIN(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) :
		im(im_, o, w, h, c, wp, pp),
		dw(static_cast<V>(w) - static_cast<V>(1.0)),
		dh(static_cast<V>(h) - static_cast<V>(1.0)) {}
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= dw && Y <= dh) {
            // Compute integer coordinates and offsets
            int y = static_cast<int>(Y);
            int x = static_cast<int>(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == dw) {
                x = x > 0 ? --x : 0;
                u = static_cast<V>(1.0);
            }
            if (Y == dh) {
                y = y > 0 ? --y : 0;
                v = static_cast<V>(1.0);
            }
            // Compute the linear index
            int k = x + y * im.WidthPitch();
            // For each image channel...
            for (int c = 0; c < im.Channels(); ++c, B += out_pitch, k += im.PlanePitch()) {
                // Do the interpolation
                V out = im.at(k) + (im.at(k+1) - im.at(k)) * u;
                out += ((im.at(k+im.WidthPitch()) - out) + (im.at(k+im.WidthPitch()+1) - im.at(k+im.WidthPitch())) * u) * v;
                *B = saturate_cast<U, V>(out);
            }
		} else {
			// Out of bounds
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch)
				*B = im.OOBV();
		}
    }
    
    // Lookup value and gradient function
    inline void lookup_grad(U *B, U *G, const V X, const V Y, const int out_pitch=1, const int grad_pitch=2) {
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= dw && Y <= dh) {
            // Compute integer coordinates and offsets
            int y = static_cast<int>(Y);
            int x = static_cast<int>(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == dw) {
                x = x > 0 ? --x : 0;
                u = static_cast<V>(1.0);
            }
            if (Y == dh) {
                y = y > 0 ? --y : 0;
                v = static_cast<V>(1.0);
            }
            // Compute the linear index
            int k = x + y * im.WidthPitch();
            // For each image channel...
            for (int c = 0; c < im.Channels(); ++c, B += out_pitch, G += grad_pitch, k += im.PlanePitch()) {
                // Sample the image
                V v00 = im.at(k);
                V v10 = im.at(k+1);
                V v01 = im.at(k+im.WidthPitch());
                V v11 = im.at(k+im.WidthPitch()+1);
                // Linearly interpolate
                V d0 = v10 - v00;
                V d1 = v11 - v01;
                V out = v00 + d0 * u;
                out += ((v01 - out) + d1 * u) * v;
                *B = saturate_cast<U, V>(out);
                d1 -= d0;
                G[1] = saturate_cast<U, V>(d0 + d1 * v);
                G[0] = saturate_cast<U, V>(v01 - v00 + d1 * u);
            }
		} else {
			// Out of bounds
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch, G += grad_pitch) {
                *B = im.OOBV();
                G[0] = static_cast<U>(0);
                G[1] = static_cast<U>(0);
            }
		}
    }

private:
	const IM_BASE<T, U, V> im;
	const V dw;
	const V dh;
};

// This shifted interpolation method is described in:
// "Linear Interpolation Revitalized", T. Blu, P. Thevenaz, M. Unser.
// IEEE Transactions on Image Processing, vol. 13, no. 5, May 2004.
#define __SHIFT_FACTOR__ 0.21132486540518713 // 0.5 * (1 - sqrt(3) / 3)
template<class T, class U, class V> class IM_SHIFT
{
public:
    // Constructor
	IM_SHIFT(const T *im_, U o, int w, int h, int c = 1, int wp = 0, int pp = 0) :
		storage(w*h*c),
		im(storage.data(), o, w, h, c, 0, 0),
		dw(static_cast<V>(w) - static_cast<V>(1.0)),
		dh(static_cast<V>(h) - static_cast<V>(1.0))
	{
		if (wp == 0)
			wp = w;
		if (pp == 0)
			pp = wp * h;
        
		// Convert the image
        const V factor = static_cast<V>(1.0 / (1.0 - __SHIFT_FACTOR__));
        const V pole = static_cast<V>(__SHIFT_FACTOR__ / (__SHIFT_FACTOR__ - 1.0));
        
        // Conversion along rows
        for (c = 0; c < im.Channels(); ++c) {
            for (h = 0; h < im.Height(); ++h) {
                V *im_ptr = &storage[c*im.PlanePitch()+h*im.WidthPitch()];
                const T *im__ptr = &im_[c*pp+h*wp];
                im_ptr[0] = static_cast<V>(im__ptr[0]);
                for (w = 1; w < im.Width(); ++w)
                    im_ptr[w] = pole * im_ptr[w-1] + factor * static_cast<V>(im__ptr[w]);
            }
        }
        
        // In place conversion along columns
        for (c = 0; c < im.Channels(); ++c) {
            for (w = 0; w < im.Width(); ++w) {
                V *im_ptr = &storage[c*im.PlanePitch()+w];
                for (h = 1; h < im.Height(); ++h, im_ptr += im.WidthPitch())
                    im_ptr[im.WidthPitch()] = pole * im_ptr[0] + factor * im_ptr[im.WidthPitch()];
            }
        }
    }
    
    // Lookup function
    inline void lookup(U *B, V X, V Y, const int out_pitch=1) {
        // Do the shift
        X -= static_cast<V>(__SHIFT_FACTOR__);
        Y -= static_cast<V>(__SHIFT_FACTOR__);
        // Check in bounds
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= dw && Y <= dh) {
            // Compute integer coordinates and offsets
            int y = static_cast<int>(Y);
            int x = static_cast<int>(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == dw) {
                x = x > 0 ? --x : 0;
                u = static_cast<V>(1.0);
            }
            if (Y == dh) {
                y = y > 0 ? --y : 0;
                v = static_cast<V>(1.0);
            }
            // Compute the linear index
            int k = x + y * im.WidthPitch();
            // For each image channel...
            for (int c = 0; c < im.Channels(); ++c, B += out_pitch, k += im.PlanePitch()) {
                // Do the interpolation
                V out = im.at(k) + (im.at(k+1) - im.at(k)) * u;
                out += ((im.at(k+im.WidthPitch()) - out) + (im.at(k+im.WidthPitch()+1) - im.at(k+im.WidthPitch())) * u) * v;
                *B = saturate_cast<U, V>(out);
            }
		} else {
			// Out of bounds
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch)
				*B = im.OOBV();
		}
    }
    
    // Lookup value and gradient function
    inline void lookup_grad(U *B, U *G, V X, V Y, const int out_pitch=1, const int grad_pitch=2) {
        // Do the shift
        X -= static_cast<V>(__SHIFT_FACTOR__);
        Y -= static_cast<V>(__SHIFT_FACTOR__);
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= dw && Y <= dh) {
            // Compute integer coordinates and offsets
            int y = static_cast<int>(Y);
            int x = static_cast<int>(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == dw) {
                x = x > 0 ? --x : 0;
                u = static_cast<V>(1.0);
            }
            if (Y == dh) {
                y = y > 0 ? --y : 0;
                v = static_cast<V>(1.0);
            }
            // Compute the linear index
            int k = x + y * im.WidthPitch();
            // For each image channel...
            for (int c = 0; c < im.Channels(); ++c, B += out_pitch, G += grad_pitch, k += im.PlanePitch()) {
                // Sample the image
                V v00 = im.at(k);
                V v10 = im.at(k+1);
                V v01 = im.at(k+im.WidthPitch());
                V v11 = im.at(k+im.WidthPitch()+1);
                // Linearly interpolate
                V d0 = v10 - v00;
                V d1 = v11 - v01;
                V out = v00 + d0 * u;
                out += ((v01 - out) + d1 * u) * v;
                *B = saturate_cast<U, V>(out);
                d1 -= d0;
                G[1] = saturate_cast<U, V>(d0 + d1 * v);
                G[0] = saturate_cast<U, V>(v01 - v00 + d1 * u);
            }
		} else {
			// Out of bounds
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch, G += grad_pitch) {
                *B = im.OOBV();
                G[0] = static_cast<U>(0);
                G[1] = static_cast<U>(0);
            }
		}
    }
    
    inline int Channels() const { return im.Channels(); }
    inline void SetOOBV(U oobv) { im.SetOOBV(oobv); }

private:
	std::vector<V> storage;
	IM_BASE<V, U, V> im;
	const V dw;
	const V dh;
};

// N TAP FILTERS
#define __PI__ 3.14159265358979323846264338327950288
template <int N>
struct lanczos { // N tap - described here: https://en.wikipedia.org/wiki/Lanczos_resampling
    static constexpr int a = N / 2;
    template <typename V> V               operator()(V x) { if (x == static_cast<V>(0)) return static_cast<V>(1.0); x *= static_cast<V>(__PI__); return sin(x) * sin(x / a) * a / (x*x); }
    template <typename V> std::pair<V, V> operator[](V x) { // Derivative
        if (x == static_cast<V>(0))
            return std::make_pair<V, V>(static_cast<V>(1.0), static_cast<V>(0.0));
        V px = x * static_cast<V>(__PI__);
        V spx = sin(px);
        V cpx = cos(px);
        V pxa = x / a;
        V spxa = sin(pxa);
        V cpxa = cos(pxa);
        return std::make_pair<V, V>(spx * spxa * a / (px * px), (spx * cpxa + a * cpx * spxa - 2 * a * spx * spxa / px) / (px * x)); 
    }
};
struct magic { // Magic kernel (3 tap) - described here: http://www.johncostella.com/magic/
    template <typename V> V               operator()(V x) { return (std::abs(x) <= static_cast<V>(0.5)) ?                     (static_cast<V>(0.75) - (x * x))                           :                     (static_cast<V>(0.5) * (x * x - static_cast<V>(3) * std::abs(x) + static_cast<V>(2.25))); }
    template <typename V> std::pair<V, V> operator[](V x) { return (std::abs(x) <= static_cast<V>(0.5)) ? std::make_pair<V, V>(static_cast<V>(0.75) - (x * x), static_cast<V>(-2.0) * x) : std::make_pair<V, V>(static_cast<V>(0.5) * (x * x - static_cast<V>(3) * std::abs(x) + static_cast<V>(2.25)), x + (x < static_cast<V>(0) ? static_cast<V>(1.5) : -static_cast<V>(1.5))); } // Derivative
};

template<class T, class U, class V, int N, typename filter> class IM_NTAP
{
public:
    static constexpr int offset = (N - 1) / 2;
    // Constructor
    IM_NTAP(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) :
		im(im_, o, w, h, c, wp, pp),
		dw(static_cast<V>(w) - static_cast<V>(1.0)),
		dh(static_cast<V>(h) - static_cast<V>(1.0)) {}
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= static_cast<V>(0.0) && X <= dw && Y >= static_cast<V>(0.0) && Y <= dh) {
            // N tap interpolation
            // Compute the filter values and lookup indices
            V xf[N], yf[N], x_val, y_val;
            int ind[N][N];
			int x = static_cast<int>(X);
			int y = static_cast<int>(Y);
            int c, d;
            filter f;
            for (c = 0; c < N; ++c) {
                xf[c] = f(X - static_cast<V>(x + c - offset));
                yf[c] = f(Y - static_cast<V>(y + c - offset));
                for (d = 0; d < N; ++d)
                    ind[c][d] = im.ind_symmetric(x + d - offset, y + c - offset);
            }
            // Sample the image
			for (c = 0; c < im.Channels(); ++c, B += out_pitch) {
				y_val = static_cast<V>(0);
                for (int y = 0; y < N; ++y) {
                    x_val = static_cast<V>(0);
                    for (int x = 0; x < N; ++x)
                        x_val += im.at(ind[y][x] + c * im.PlanePitch()) * xf[x];
					y_val += yf[y] * x_val;
                }
				*B = saturate_cast<U, V>(y_val);
			}
		} else {
			// Out of bounds
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch)
				*B = im.OOBV();
		}
    }
    
    // Lookup value and gradient function
    inline void lookup_grad(U *B, U *G, const V X, const V Y, const int out_pitch=1, const int grad_pitch=2) {
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= dw && Y <= dh) {
            // N tap interpolation
            // Compute the filter values and lookup indices
            std::pair<V, V> xf[N], yf[N];
            V x_val, y_val, im_val, x_deriv_, x_deriv, y_deriv;
            int ind[N][N];
			int x = static_cast<int>(X);
			int y = static_cast<int>(Y);
            int c, d;
            filter f;
            for (c = 0; c < N; ++c) {
                xf[c] = f[X - static_cast<V>(x + c - offset)];
                yf[c] = f[Y - static_cast<V>(y + c - offset)];
                for (d = 0; d < N; ++d)
                    ind[c][d] = im.ind_symmetric(x + d - offset, y + c - offset);
            }
            // Sample the image
			for (c = 0; c < im.Channels(); ++c, B += out_pitch, G += grad_pitch) {
				y_val = static_cast<V>(0);
                x_deriv = static_cast<V>(0);
                y_deriv = static_cast<V>(0);
                for (int y = 0; y < N; ++y) {
                    x_val = static_cast<V>(0);
                    x_deriv_ = static_cast<V>(0);
                    for (int x = 0; x < N; ++x) {
                        im_val = im.at(ind[y][x] + c * im.PlanePitch());
                        x_val += im_val * xf[x].first;
                        x_deriv_ += im_val * xf[x].second;
                    }
					y_val += x_val * yf[y].first;
					x_deriv += x_deriv_ * yf[y].first;
					y_deriv += x_val * yf[y].second;
                }
				*B = saturate_cast<U, V>(y_val);
                G[0] = saturate_cast<U, V>(x_deriv);
                G[1] = saturate_cast<U, V>(y_deriv);
			}
		} else {
			// Out of bounds
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch, G += grad_pitch) {
                *B = im.OOBV();
                G[0] = static_cast<U>(0);
                G[1] = static_cast<U>(0);
            }
		}
    }

private:
	const IM_BASE<T, U, V> im;
	const V dw;
	const V dh;
};

template<class T, class U, class V> class IM_CUB
{
public:
    // Constructor
    IM_CUB(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) :
		im(im_, o, w, h, c, wp, pp),
		dw(static_cast<V>(w) - static_cast<V>(2.0)),
		dh(static_cast<V>(h) - static_cast<V>(2.0)) {}
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= static_cast<V>(2.0) && X < dw && Y >= static_cast<V>(2.0) && Y < dh) {
            // Bicubicly interpolate
            V b[4], d[4], u[3], v[3];
			int x = static_cast<int>(X);
			int y = static_cast<int>(Y);
			u[0] = X - x;
			v[0] = Y - y;
			u[1] = u[0] * u[0];
			v[1] = v[0] * v[0];
			u[2] = u[1] * u[0];
			v[2] = v[1] * v[0];
			int k = x - 1 + (y - 1) * im.WidthPitch();
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch, k += im.PlanePitch()) {
				V a;
                for (int m = 0, n = k; m < 4; ++m, n += im.WidthPitch()) {
					d[0] = im.at(n+0);
					d[1] = im.at(n+1);
					d[2] = im.at(n+2);
					d[3] = im.at(n+3);
					a = (d[3] + d[1]) - (d[2] + d[0]);
					b[m] = u[2] * a + u[1] * ((d[0] - d[1]) - a) + u[0] * (d[2] - d[0]) + d[1];
				}
				a = (b[3] + b[1]) - (b[2] + b[0]);
				*B = saturate_cast<U, V>(v[2] * a + v[1] * ((b[0] - b[1]) - a) + v[0] * (b[2] - b[0]) + b[1]);
			}
		} else {
			// Out of bounds
			for (int c = 0; c < im.Channels(); ++c, B += out_pitch)
				*B = im.OOBV();
		}
    }
    
    inline void lookup_grad(U *B, U *G, const V X, const V Y, const int out_pitch=1, const int grad_pitch=2) {}

private:
	const IM_BASE<T, U, V> im;
	const V dw;
	const V dh;
};

#endif
