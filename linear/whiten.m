%WHITEN Transform data to be zero mean and identity covariance
%
%   [Y, T] = whiten(X, [epsilon])
%
%IN:
%   X - MxN array of N vectors od dimension M to be whitened
%   epsilon - scalar value to add to eigen values to avoid amplifying
%             noise. Default: 1e-4.
%
%   Y - MxN array of whitened data.
%   T - (M+1)x(M+1) matrix for the transformation [Y; 1] = T * [X; 1].

function [X, T] = whiten(X, epsilon)
% Set the default epsilon
if nargin < 2
    epsilon = 1e-4;
end

% Subtract the mean
mu = mean(X, 2);
X = bsxfun(@minus, X, mu);

% Whitening the covariance
[V, D] = svd(X * X');
T = V' * diag(sqrt((size(X, 2) - 1) ./ (diag(D) + epsilon))) * V;
X = T * X;

% Construct the transform
mu = T * -mu;
T(end+1,end+1) = 1;
T(1:end-1,end) = mu;
end
