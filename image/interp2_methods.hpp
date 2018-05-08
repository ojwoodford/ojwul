#ifndef __INTERP2_METHODS_H__
#define __INTERP2_METHODS_H__
#include <limits> // For std::numeric_limits
#include <utility> // For std::pair
#include <cstdlib> // For std::abs

// These methods all use the following conventions:
//    1. The centre of the top left pixel is at 0,0
//    2. The image is stored in memory row-major.

template<class T, class U, class V> class IM_BASE
{
public:
    // Image dimensions
	inline int Width() const { return width; }
	inline int Height() const { return height; }
	inline int Channels() const { return nchannels; }
    
    // Other helpers
	inline void SetOOBV(U new_oobv) { oobv = new_oobv; }
    
protected:    
    const T* im;
    V* im_copy;
    int height;
    int width;
    int nchannels;
    int width_pitch;
    int plane_pitch;
    U oobv;
    V dw;
    V dh;
    
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
};

// Function for correct rounding
template <typename U, typename T>
static inline U saturate_cast(T val)
{
    U out = static_cast<U>(val);
    if (std::numeric_limits<U>::is_integer && !std::numeric_limits<T>::is_integer) {
        if (std::numeric_limits<U>::is_signed)
            return out > static_cast<U>(0) ? (out > std::numeric_limits<U>::max() ? std::numeric_limits<U>::max() : out + static_cast<U>(0.5)) : (out < std::numeric_limits<U>::min() ? std::numeric_limits<U>::min() : out - static_cast<U>(0.5));
        else
            return out > static_cast<U>(0) ? (out > std::numeric_limits<U>::max() ? std::numeric_limits<U>::max() : out + static_cast<U>(0.5)) : static_cast<U>(0);
    }
    return out;
}

template<class T, class U, class V> class IM_NEAR : public IM_BASE<T, U, V>
{
public:
    // Constructor
    IM_NEAR(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) {
        IM_BASE<T, U, V>::im = im_;
        IM_BASE<T, U, V>::oobv = o;
        IM_BASE<T, U, V>::width = w;
        IM_BASE<T, U, V>::dw = static_cast<V>(w) - static_cast<V>(0.5);
        IM_BASE<T, U, V>::height = h;
        IM_BASE<T, U, V>::dh = static_cast<V>(h) - static_cast<V>(0.5);
        IM_BASE<T, U, V>::nchannels = c;
        IM_BASE<T, U, V>::width_pitch = wp == 0 ? IM_BASE<T, U, V>::width : wp;
        IM_BASE<T, U, V>::plane_pitch = pp == 0 ? IM_BASE<T, U, V>::width_pitch * IM_BASE<T, U, V>::height : pp;
    }
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= -0.5 && X < IM_BASE<T, U, V>::dw && Y >= -0.5 && Y < IM_BASE<T, U, V>::dh) {
			// Find nearest neighbour
			int k = static_cast<int>(X + static_cast<V>(0.5)) + static_cast<int>(Y + static_cast<V>(0.5)) * IM_BASE<T, U, V>::width_pitch;
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, k += IM_BASE<T, U, V>::plane_pitch)
				*B = saturate_cast<U, T>(IM_BASE<T, U, V>::im[k]);
		} else {
			// Out of bounds
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch)
				*B = IM_BASE<T, U, V>::oobv;
		}
    }
    
    // Lookup value and gradient function
    inline void lookup_grad(U *B, U *G, const V X, const V Y, const int out_pitch=1, const int grad_pitch=2) {
        if (X >= -0.5 && X < IM_BASE<T, U, V>::dw && Y >= -0.5 && Y < IM_BASE<T, U, V>::dh) {
			// Find nearest neighbour
            int y = static_cast<int>(Y + static_cast<V>(0.5));
            int x = static_cast<int>(X + static_cast<V>(0.5));
			int k = x + y * IM_BASE<T, U, V>::width_pitch;
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, G += grad_pitch, k += IM_BASE<T, U, V>::plane_pitch)
            {
				*B = saturate_cast<U, T>(IM_BASE<T, U, V>::im[k]);
                if (y == 0)
                    G[0] = saturate_cast<U, V>(static_cast<V>(IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch]) - static_cast<V>(IM_BASE<T, U, V>::im[k]));
                else if (y == IM_BASE<T, U, V>::height-1)
                    G[0] = saturate_cast<U, V>(static_cast<V>(IM_BASE<T, U, V>::im[k]) - static_cast<V>(IM_BASE<T, U, V>::im[k-IM_BASE<T, U, V>::width_pitch]));
                else
                    G[0] = saturate_cast<U, V>(static_cast<V>(0.5) * (static_cast<V>(IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch]) - static_cast<V>(IM_BASE<T, U, V>::im[k-IM_BASE<T, U, V>::width_pitch])));
                if (x == 0)
                    G[1] = saturate_cast<U, V>(static_cast<V>(IM_BASE<T, U, V>::im[k+1]) - static_cast<V>(IM_BASE<T, U, V>::im[k]));
                else if (x == IM_BASE<T, U, V>::width-1)
                    G[1] = saturate_cast<U, V>(static_cast<V>(IM_BASE<T, U, V>::im[k]) - static_cast<V>(IM_BASE<T, U, V>::im[k-1]));
                else
                    G[1] = saturate_cast<U, V>(static_cast<V>(0.5) * (static_cast<V>(IM_BASE<T, U, V>::im[k+1]) - static_cast<V>(IM_BASE<T, U, V>::im[k-1])));
            }
		} else {
			// Out of bounds
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, G += grad_pitch) {
                *B = IM_BASE<T, U, V>::oobv;
                G[0] = static_cast<U>(0);
                G[1] = static_cast<U>(0);
            }
		}
    }
};

template<class T, class U, class V> class IM_LIN : public IM_BASE<T, U, V>
{
public:
    // Constructor
    IM_LIN(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) {
        IM_BASE<T, U, V>::im = im_;
        IM_BASE<T, U, V>::oobv = o;
        IM_BASE<T, U, V>::width = w;
        IM_BASE<T, U, V>::dw = static_cast<V>(w) - static_cast<V>(1.0);
        IM_BASE<T, U, V>::height = h;
        IM_BASE<T, U, V>::dh = static_cast<V>(h) - static_cast<V>(1.0);
        IM_BASE<T, U, V>::nchannels = c;
        IM_BASE<T, U, V>::width_pitch = wp == 0 ? IM_BASE<T, U, V>::width : wp;
        IM_BASE<T, U, V>::plane_pitch = pp == 0 ? IM_BASE<T, U, V>::width_pitch * IM_BASE<T, U, V>::height : pp;
    }
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= IM_BASE<T, U, V>::dw && Y <= IM_BASE<T, U, V>::dh) {
            // Compute integer coordinates and offsets
            int y = static_cast<int>(Y);
            int x = static_cast<int>(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == IM_BASE<T, U, V>::dw) {
                x = x > 0 ? --x : 0;
                u = static_cast<V>(1.0);
            }
            if (Y == IM_BASE<T, U, V>::dh) {
                y = y > 0 ? --y : 0;
                v = static_cast<V>(1.0);
            }
            // Compute the linear index
            int k = x + y * IM_BASE<T, U, V>::width_pitch;
            // For each image channel...
            for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, k += IM_BASE<T, U, V>::plane_pitch) {
                // Do the interpolation
                V out = static_cast<V>(IM_BASE<T, U, V>::im[k]) + (static_cast<V>(IM_BASE<T, U, V>::im[k+1]) - static_cast<V>(IM_BASE<T, U, V>::im[k])) * u;
                out += ((static_cast<V>(IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch]) - out) + (static_cast<V>(IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch+1]) - static_cast<V>(IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch])) * u) * v;
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
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= IM_BASE<T, U, V>::dw && Y <= IM_BASE<T, U, V>::dh) {
            // Compute integer coordinates and offsets
            int y = static_cast<int>(Y);
            int x = static_cast<int>(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == IM_BASE<T, U, V>::dw) {
                x = x > 0 ? --x : 0;
                u = static_cast<V>(1.0);
            }
            if (Y == IM_BASE<T, U, V>::dh) {
                y = y > 0 ? --y : 0;
                v = static_cast<V>(1.0);
            }
            // Compute the linear index
            int k = x + y * IM_BASE<T, U, V>::width_pitch;
            // For each image channel...
            for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, G += grad_pitch, k += IM_BASE<T, U, V>::plane_pitch) {
                // Sample the image
                V v00 = static_cast<V>(IM_BASE<T, U, V>::im[k]);
                V v10 = static_cast<V>(IM_BASE<T, U, V>::im[k+1]);
                V v01 = static_cast<V>(IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch]);
                V v11 = static_cast<V>(IM_BASE<T, U, V>::im[k+IM_BASE<T, U, V>::width_pitch+1]);
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
                G[0] = static_cast<U>(0);
                G[1] = static_cast<U>(0);
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
        IM_BASE<T, U, V>::dw = static_cast<V>(w) - static_cast<V>(1.0);
        IM_BASE<T, U, V>::height = h;
        IM_BASE<T, U, V>::dh = static_cast<V>(h) - static_cast<V>(1.0);
        IM_BASE<T, U, V>::nchannels = c;
        wp = wp == 0 ? IM_BASE<T, U, V>::width : wp;
        pp = pp == 0 ? wp * IM_BASE<T, U, V>::height : pp;
        IM_BASE<T, U, V>::width_pitch = IM_BASE<T, U, V>::width;
        IM_BASE<T, U, V>::plane_pitch = IM_BASE<T, U, V>::width_pitch * IM_BASE<T, U, V>::height;
        
        // Allocate memory for the new coefficients
        IM_BASE<T, U, V>::im_copy = new V[IM_BASE<T, U, V>::plane_pitch*IM_BASE<T, U, V>::nchannels];
        
        const V factor = static_cast<V>(1.0 / (1.0 - __SHIFT_FACTOR__));
        const V pole = static_cast<V>(__SHIFT_FACTOR__ / (__SHIFT_FACTOR__ - 1.0));
        
        // Conversion along rows
        for (c = 0; c < IM_BASE<T, U, V>::nchannels; ++c) {
            for (h = 0; h < IM_BASE<T, U, V>::height; ++h) {
                V *im_ptr = &IM_BASE<T, U, V>::im_copy[c*IM_BASE<T, U, V>::plane_pitch+h*IM_BASE<T, U, V>::width_pitch];
                const T *im__ptr = &im_[c*pp+h*wp];
                im_ptr[0] = static_cast<V>(im__ptr[0]);
                for (w = 1; w < IM_BASE<T, U, V>::width; ++w)
                    im_ptr[w] = pole * im_ptr[w-1] + factor * static_cast<V>(im__ptr[w]);
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
        X -= static_cast<V>(__SHIFT_FACTOR__);
        Y -= static_cast<V>(__SHIFT_FACTOR__);
        // Check in bounds
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= IM_BASE<T, U, V>::dw && Y <= IM_BASE<T, U, V>::dh) {
            // Compute integer coordinates and offsets
            int y = static_cast<int>(Y);
            int x = static_cast<int>(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == IM_BASE<T, U, V>::dw) {
                x = x > 0 ? --x : 0;
                u = static_cast<V>(1.0);
            }
            if (Y == IM_BASE<T, U, V>::dh) {
                y = y > 0 ? --y : 0;
                v = static_cast<V>(1.0);
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
        X -= static_cast<V>(__SHIFT_FACTOR__);
        Y -= static_cast<V>(__SHIFT_FACTOR__);
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= IM_BASE<T, U, V>::dw && Y <= IM_BASE<T, U, V>::dh) {
            // Compute integer coordinates and offsets
            int y = static_cast<int>(Y);
            int x = static_cast<int>(X);
            V v = Y - y;
            V u = X - x;
            // Check for boundary cases
            if (X == IM_BASE<T, U, V>::dw) {
                x = x > 0 ? --x : 0;
                u = static_cast<V>(1.0);
            }
            if (Y == IM_BASE<T, U, V>::dh) {
                y = y > 0 ? --y : 0;
                v = static_cast<V>(1.0);
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
                G[0] = static_cast<U>(0);
                G[1] = static_cast<U>(0);
            }
		}
    }
};

// N TAP FILTERS
#define __PI__ 3.14159265358979323846264338327950288
template <int N>
struct lanczos { // N tap - described here: https://en.wikipedia.org/wiki/Lanczos_resampling
    template <typename V> V               operator()(V x) { if (x == static_cast<V>(0)) return static_cast<V>(1.0); x *= static_cast<V>(__PI__); return sin(x) * sin(x / (N-1)) * (N-1) / (x*x); }
    template <typename V> std::pair<V, V> operator[](V x) { // Derivative
        if (x == static_cast<V>(0))
            return std::make_pair<V, V>(static_cast<V>(1.0), static_cast<V>(0.0));
        V px = x * static_cast<V>(__PI__);
        V spx = sin(px);
        V cpx = cos(px);
        V pxa = x / (N - 1);
        V spxa = sin(pxa);
        V cpxa = cos(pxa);
        return std::make_pair<V, V>(spx * spxa * (N-1) / (px * px), (spx * cpxa + (N - 1) * cpx * spxa - 2 * (N - 1) * spx * spxa / px) / (px * x)); 
    }
};
struct magic { // Magic kernel (3 tap) - described here: http://johncostella.webs.com/magic/
    template <typename V> V               operator()(V x) { return (std::abs(x) <= static_cast<V>(0.5)) ?                     (static_cast<V>(0.75) - (x * x))                           :                     (static_cast<V>(0.5) * (x * x - static_cast<V>(3) * std::abs(x) + static_cast<V>(2.25))); }
    template <typename V> std::pair<V, V> operator[](V x) { return (std::abs(x) <= static_cast<V>(0.5)) ? std::make_pair<V, V>(static_cast<V>(0.75) - (x * x), static_cast<V>(-2.0) * x) : std::make_pair<V, V>(static_cast<V>(0.5) * (x * x - static_cast<V>(3) * std::abs(x) + static_cast<V>(2.25)), x + (x < static_cast<V>(0) ? static_cast<V>(1.5) : -static_cast<V>(1.5))); } // Derivative
};

template<class T, class U, class V, int N, typename filter> class IM_NTAP : public IM_BASE<T, U, V>
{
public:
    // Constructor
    IM_NTAP(const T *im_, U o, int w, int h, int c=1, int wp=0, int pp=0) {
        IM_BASE<T, U, V>::im = im_;
        IM_BASE<T, U, V>::oobv = o;
        IM_BASE<T, U, V>::width = w;
        IM_BASE<T, U, V>::dw = static_cast<V>(w) - static_cast<V>(1.0);
        IM_BASE<T, U, V>::height = h;
        IM_BASE<T, U, V>::dh = static_cast<V>(h) - static_cast<V>(1.0);
        IM_BASE<T, U, V>::nchannels = c;
        IM_BASE<T, U, V>::width_pitch = wp == 0 ? IM_BASE<T, U, V>::width : wp;
        IM_BASE<T, U, V>::plane_pitch = pp == 0 ? IM_BASE<T, U, V>::width_pitch * IM_BASE<T, U, V>::height : pp;
    }
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= static_cast<V>(0.0) && X <= IM_BASE<T, U, V>::dw && Y >= static_cast<V>(0.0) && Y <= IM_BASE<T, U, V>::dh) {
            // N tap interpolation
            // Compute the filter values and lookup indices
            V xf[N], yf[N], x_val, y_val;
            const T* im_;
            int ind[N][N];
			int x = static_cast<int>(X);
			int y = static_cast<int>(Y);
            int c, d;
            filter f;
            for (c = 0; c < N; ++c) {
                xf[c] = f(X - static_cast<V>(x + c - (N - 1) / 2));
                yf[c] = f(Y - static_cast<V>(y + c - (N - 1) / 2));
                for (d = 0; d < N; ++d)
                    ind[c][d] = IM_BASE<T, U, V>::ind_symmetric(x + d - (N - 1) / 2, y + c - (N - 1) / 2);
            }
            // Sample the image
            im_ = IM_BASE<T, U, V>::im;
			for (c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, im_ += IM_BASE<T, U, V>::plane_pitch) {
				y_val = static_cast<V>(0);
                for (int y = 0; y < N; ++y) {
                    x_val = static_cast<V>(0);
                    for (int x = 0; x < N; ++x)
                        x_val += static_cast<V>(im_[ind[y][x]] * xf[x]);
					y_val += yf[y] * x_val;
                }
				*B = saturate_cast<U, V>(y_val);
			}
		} else {
			// Out of bounds
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch)
				*B = IM_BASE<T, U, V>::oobv;
		}
    }
    
    // Lookup value and gradient function
    inline void lookup_grad(U *B, U *G, const V X, const V Y, const int out_pitch=1, const int grad_pitch=2) {
        if (X >= static_cast<V>(0.0) && Y >= static_cast<V>(0.0) && X <= IM_BASE<T, U, V>::dw && Y <= IM_BASE<T, U, V>::dh) {
            // N tap interpolation
            // Compute the filter values and lookup indices
            std::pair<V, V> xf[N], yf[N];
            V x_val, y_val, im_val, x_deriv_, x_deriv, y_deriv;
            const T* im_;
            int ind[N][N];
			int x = static_cast<int>(X);
			int y = static_cast<int>(Y);
            int c, d;
            filter f;
            for (c = 0; c < N; ++c) {
                xf[c] = f[X - static_cast<V>(x + c - (N - 1) / 2)];
                yf[c] = f[Y - static_cast<V>(y + c - (N - 1) / 2)];
                for (d = 0; d < N; ++d)
                    ind[c][d] = IM_BASE<T, U, V>::ind_symmetric(x + d - (N - 1) / 2, y + c - (N - 1) / 2);
            }
            // Sample the image
            im_ = IM_BASE<T, U, V>::im;
			for (c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, G += grad_pitch, im_ += IM_BASE<T, U, V>::plane_pitch) {
				y_val = static_cast<V>(0);
                x_deriv = static_cast<V>(0);
                y_deriv = static_cast<V>(0);
                for (int y = 0; y < N; ++y) {
                    x_val = static_cast<V>(0);
                    x_deriv_ = static_cast<V>(0);
                    for (int x = 0; x < N; ++x) {
                        im_val = static_cast<V>(im_[ind[y][x]]);
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
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, G += grad_pitch) {
                *B = IM_BASE<T, U, V>::oobv;
                G[0] = static_cast<U>(0);
                G[1] = static_cast<U>(0);
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
        IM_BASE<T, U, V>::dw = static_cast<V>(w) - static_cast<V>(2.0);
        IM_BASE<T, U, V>::height = h;
        IM_BASE<T, U, V>::dh = static_cast<V>(h) - static_cast<V>(2.0);
        IM_BASE<T, U, V>::nchannels = c;
        IM_BASE<T, U, V>::width_pitch = wp == 0 ? IM_BASE<T, U, V>::width : wp;
        IM_BASE<T, U, V>::plane_pitch = pp == 0 ? IM_BASE<T, U, V>::width_pitch * IM_BASE<T, U, V>::height : pp;
    }
    
    // Lookup function
    inline void lookup(U *B, const V X, const V Y, const int out_pitch=1) {
        if (X >= static_cast<V>(2.0) && X < IM_BASE<T, U, V>::dw && Y >= static_cast<V>(2.0) && Y < IM_BASE<T, U, V>::dh) {
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
			int k = x - 1 + (y - 1) * IM_BASE<T, U, V>::width_pitch;
			for (int c = 0; c < IM_BASE<T, U, V>::nchannels; ++c, B += out_pitch, k += IM_BASE<T, U, V>::plane_pitch) {
				V a;
                for (int m = 0, n = k; m < 4; ++m, n += IM_BASE<T, U, V>::width_pitch) {
					d[0] = static_cast<V>(IM_BASE<T, U, V>::im[n+0]);
					d[1] = static_cast<V>(IM_BASE<T, U, V>::im[n+1]);
					d[2] = static_cast<V>(IM_BASE<T, U, V>::im[n+2]);
					d[3] = static_cast<V>(IM_BASE<T, U, V>::im[n+3]);
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
