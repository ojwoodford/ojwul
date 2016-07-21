%PCA Principal component analysis
%
%   [T, Y, V] = pca(X)
%
% IN:
%    X - MxN matrix of N vectors of dimension M.
%
% OUT:
%    T - Mx(M+1) Projection matrix for PCA transformation
%    Y - MxN matrix of transformed vectors, where Y = T * [X; ones(1, N].
%    V - Mx1 list of eigen values associated with each dimension.

function [T, X, V] = pca(X)

% Compute the mean and subtract
t = mean(X, 2);
X = bsxfun(@minus, X, t);

% Do the economy svd, to avoid memory issues
[R, V] = svd(X, 'econ');
R = R';

% Compose the transform
T = [R -(R * t)];

if nargout > 1
    % Apply the transform
    X = R * X;
end

if nargout > 2
    % Square the singular values to give eigenvalues
    V = diag(V);
    V = V .* V;
end
end