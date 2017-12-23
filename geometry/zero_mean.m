%ZERO_MEAN Subtract the means from a set of vectors to make them zero mean
%
%   Y = zero_mean(X)
%   Y = zero_mean(X, dim)
%
% Subtract the mean along a specified dimension from all vectors in an 
% array, making them zero mean.
%
%IN:
%   X - Array containing vectors to zero mean.
%   dim - Dimension along which to zero-mean X. Default: first
%         non-singleton dimension of X.
%
%OUT:
%   Y - Zero-meaned X.

function X = zero_mean(X, varargin)
X = bsxfun(@minus, X, mean(X, varargin{:}));
end
