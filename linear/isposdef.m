%ISPOSDEF Check if a square matrix is positive definite
%
%   tf = isposdef(A)
%
%IN:
%   A - NxN matrix.
%
%OUT:
%   tf - boolean indicating whether A is positive definite.

function tf = isposdef(A)
[~, tf] = chol(A);
tf = (tf == 0) && (rank(A) == size(A, 1));
end