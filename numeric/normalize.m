%NORMALIZE Set vectors in an array to be of unit length
%
%   Y = normalize(X)
%   Y = normalize(X, dim)
%
% Set all the non-zero vectors in an array, along a specified dimension, to
% be of unit length.
%
%IN:
%   X - Array containing vectors to normalize.
%   dim - Dimension along which to normalize X. Default: first
%         non-singleton dimension of X.
%
%OUT:
%   Y - Normalize X.

function X = normalize(X, varargin)
normalizer = normd(X, varargin{:});
M = normalizer ~= 0;
normalizer(M) = 1 ./ normalizer(M);
X = bsxfun(@times, X, normalizer);
end
