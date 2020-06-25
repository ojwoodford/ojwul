%COL Convert an array to a column vector along a particular dimension
%
%   B = col(A, [dim])
%
%IN:
%   A - Array of any size.
%   dim - Positive integer indicating the dimension to arrange the elements
%         of A along. Default: 1.
%
%OUT:
%   B - Result of shiftdim(A(:), 1-dim).

function x = col(x, dim)
x = reshape(x, numel(x), 1);
if nargin > 1
    x = shiftdim(x, 1-dim);
end
end