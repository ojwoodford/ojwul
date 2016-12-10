%FIRST Returns indices of the first non-zero elements along the given dimension
%
%   B = first(A, [dim])
%
% For each vector along the given dimension, this function returns the
% index of the first non-zero element along that vector, or 0 if there is
% no non-zero element.
%
%IN:
%   A - MxNx... input array
%   dim - dimension of A along which to find the first non-zero element.
%         Default: first non-singleton dimension of A.
%
%OUT
%   B - 1xNx... output array (size dependent on which dimension specified).

function A = first(A, dim)
if nargin < 2
    dim = find(size(A) ~= 1, 'first', 1);
    if isempty(dim)
        dim = 1;
    end
end
if ~islogical(A)
    A = A ~= 0;
end
n = size(A, dim);
A = cumsum(A, dim);
A = sum(A > 0, dim);
A = (n + 1) - A;
A(A==n+1) = 0;
end
