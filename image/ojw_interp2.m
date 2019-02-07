%OJW_INTERP2  Fast 2d interpolation for images
%
%	V = ojw_interp2(A, X, Y)
%	V = ojw_interp2(A, X, Y, interp_mode)
%	V = ojw_interp2(A, X, Y, interp_mode, oobv)
%
% 2d interpolation on a regular grid - similar to MATLAB's interp2() but
% with much less overhead, and supports multiple channels and types. Note
% that while results of 'linear' and 'nearest' interpolation are the same
% as those of interp2(), those of cubic are not - ojw_interp2 uses a cubic
% hermite spline that is very fast to compute, unlike the natural cubic
% spline employed by interp2(), which does, however, yield a smoother
% interpolation.
%
%IN:
%	A - HxWxC double, single, uint16, int16, uint8, int8 or Logical array.
%	X - MxN horizontal offsets (1 being the centre of the first pixel).
%	Y - MxN vertical offsets (1 being the centre of the first pixel).
%	interp_mode - string, either 'nearest', 'linear', 'cubic', 'magic',
%	              '3lanczos', or '5lanczos'. Default: 'linear'.
%	oobv - 1x1 Out of bounds value. Default: NaN.
%
%OUT:
%	V - MxNxC interpolated values. Class is the same as that of oobv.
%   G - 2xMxNxC array of gradients of V in x and y directions:
%       [shiftdim(Vx, -1); shiftdim(Vy, -1)].

function varargout = ojw_interp2(varargin)
sourceList = {'ojw_interp2.cpp', '-Xopenmp'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
end