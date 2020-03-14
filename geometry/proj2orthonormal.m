%PROJ2ORTHONORMAL Project a matrix onto an orthonormal basis
%
%   B = proj2orthonormal(A)
%
%IN:
%   A - MxN matrix
%   B - Closest MxN matrix to A, where the rows are unit length and
%       orthogonal.

function A = proj2orthonormal(A)
[U, ~, V] = svd(A, 'econ');
A = U * V';
end