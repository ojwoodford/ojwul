%OJW_INTERP2_ALT  Fast 2d interpolation for images
%
%	V = ojw_interp2_alt(A, X)
%	V = ojw_interp2_alt(A, X, interp_mode)
%	V = ojw_interp2_alt(A, X, interp_mode, oobv)
%
% Wrapper to ojw_interp2 which has horizontal and vertical coordinates
% concatenated along the third dimension into one array.
%
%IN:
%	A - HxWxC double, single, uint16, int16, uint8, int8 or Logical array.
%	X - MxNx2 sample coordinates (1,1 being the centre of the top left
%	    pixel). 
%	interp_mode - string, either 'nearest', 'linear', 'cubic', 'magic',
%	              '3lanczos', or '5lanczos'. Default: 'linear'.
%	oobv - 1x1 Out of bounds value. Default: NaN.
%
%OUT:
%	V - MxNxC interpolated values. Class is the same as that of oobv.

function varargout = ojw_interp2_alt(A, X, varargin)
[varargout{1:nargout}] = ojw_interp2(A, X(:,:,1), X(:,:,2), varargin{:});
end