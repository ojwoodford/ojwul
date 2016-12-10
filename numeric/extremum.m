%EXTREMUM Compute the extreme value along a given dimension
%
%   B = extremum(A, [dim])
%
% Output the most extreme value (furthest from 0) along a specified
% dimension of an input array.
%
%IN:
%   A - Numeric input array.
%   dim - Dimension along which to compute the extremum. Default: first
%         non-singleton dimension.
%
%OUT:
%   B - Output array.

function A = extremum(A, dim)
if nargin < 2
    % Find the first non-singleton dimension
    dim = find(size(A) > 1, 1, 'first');
    if isempty(dim)
        return;
    end
end
[I, I] = max(abs(A), [], dim);
A = dimsel(A, I);