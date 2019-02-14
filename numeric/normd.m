%NORMD Compute the 2-norms of vectors in an array along a specific dimension
%
%   Y = normd(X)
%   Y = normd(X, dim)
%
% Compute the 2-norms of vectors in an array, along a specific dimension.
%
%IN:
%   X - Array containing vectors to compute the norms of.
%   dim - Dimension along which to compute the 2-norm. Default: first
%         non-singleton dimension of X.
%
%OUT:
%   Y - 2-norms of X.

function X = normd(X, varargin)
X = sqrt(sum(X .* X, varargin{:}));