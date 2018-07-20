%PLANE_FIT Weighted least squares plane fit to 3D points
%
%   N = plane_fit(X, [W])
%
% Finds the (weighted) least squares fit of a plane to 3D points, such that
% N' * homg(X) = 0.
%
%IN:
%   X - 3xM matrix of 3D points.
%   W - 1xM vector of weights per point. Default: ones(1, M).
%
%OUT:
%   N - 1x4 plane equation vector.

function N = plane_fit(X, W)
m = mean(X, 2);
X = bsxfun(@minus, X, m);
if nargin > 1
    X = bsxfun(@times, X, sqrt(W(:)'));
end
[N, ~] = svd(X, 0);
N = N(:,3)';
N(4) = -(N * m);
end
