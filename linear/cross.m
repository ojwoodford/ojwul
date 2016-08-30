%CROSS  Vector cross product.
%   C = CROSS(A,B) returns the cross product of the vectors
%   A and B.  That is, C = A x B.  A and B must be 3 element
%   vectors.
%
%   C = CROSS(A,B) returns the cross product of A and B along the
%   first dimension of length 3.
%
%   Class support for inputs A,B:
%      float: double, single
%
%   See also DOT.

function X = cross(X, Y)
X = X .* Y([2 3 1],:,:) - X([2 3 1],:,:) .* Y;
X = X([2 3 1],:,:);
end