%WHITEN_SRT Transform data to be zero mean and close to identity covariance
%
%   [Y, T] = whiten_srt(X)
%
% Apply a similarity transform to the data such that the mean is zero and
% the covariance is as close to the identity as possible.
%
%IN:
%   X - MxN array of N vectors of dimension M to be whitened
%
%   Y - (M+1)xN array of whitened data.
%   T - (M+1)x(M+1) matrix for the transformation proj(T \ homg(Y)) = X.

function [X, T] = whiten_srt(X)
% Subtract the mean
mu = mean(X, 2);
X = bsxfun(@minus, X, mu);

% Whitening the covariance
[T, s] = svd(X * X');
T = T * sqrt((size(X, 2) - 1) ./ max(diag(s)));
X = T * X;

% Construct the transform
T(:,end+1) = T * -mu;
T(3,3) = 1;
end
