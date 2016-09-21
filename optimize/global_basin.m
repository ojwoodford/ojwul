%GLOBAL_BASIN Output binary mask of watershed for the global min of array
%
%   B = global_basin(A)
%
% This function computes the binary mask of all those points in an array
% from which local minimization (e.g. gradient descent with small steps)
% would lead to the global minimum of the array. This is the watershed of
% the global minimum. This assumes that the input array is a regular
% sampling of a Euclidean cost space.
%
%IN:
%   A - N-dim array of scalar values
%
%OUT:
%   B - size(A) logical array, with true indicating this point would lead
%       to the global minimum.

function B = global_basin(A)
B = watershed(A, conndef(ndims(A), 'maximal'));
[b, b] = min(A(:));
B = B == B(b);