%ALL_FINITE Checks if all the elements in an array are finite
%
%   tf = all_finite(x)
%
%IN:
%   x - Numeric array.
%
%OUT:
%   tf - Boolean indicating whether all elements of x are finite.

function x = all_finite(x)
if issparse(x)
    [~, ~, x] = find(x);
end
x = all(isfinite(x(:)));
end