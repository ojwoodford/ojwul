function c = dot(a,b,dim)
%DOT  Vector dot product.
%   C = DOT(A,B) returns the scalar product of the vectors A and B.
%   A and B must be vectors of the same length.  When A and B are both
%   column vectors, DOT(A,B) is the same as A'*B.
%
%   DOT(A,B), for N-D arrays A and B, returns the scalar product
%   along the first non-singleton dimension of A and B. A and B must have
%   compatible sizes.
%
%   DOT(A,B,DIM) returns the scalar product of A and B in the
%   dimension DIM.
%
%   Class support for inputs A,B:
%      float: double, single
%
%   See also CROSS.

%   Copyright 1984-2011 The MathWorks, Inc. 

if isinteger(a) || isinteger(b) 
    error(message('MATLAB:dot:integerClass'));
end

% Special case: A and B are vectors and dim not supplied
if ismatrix(a) && ismatrix(b) && nargin < 3
   if min(size(a)) == 1, a = a(:); end
   if min(size(b)) == 1, b = b(:); end
end;

% Check dimensions
sza = size(a);
szb = size(b);
sza(end+1:numel(sz(b))) = 1;
szb(end+1:numel(sz(a))) = 1;
if any(sza ~= szb & sza ~= 1 & szb ~= 1)
   error(message('MATLAB:dot:InputSizeMismatch'));
end

if nargin == 2
    dim = find(sza ~= 1 & szb ~= 1, 1);
else
    if ~isnumeric(dim)
        error(message('MATLAB:getdimarg:dimensionMustBePositiveInteger'));
    end
end
c = sum(conj(a) .* b, dim);
end
