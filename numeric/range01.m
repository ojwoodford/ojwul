%RANGE01 Apply gain and bias so range of data is exactly [0, 1]
%
%   Y = range01(X, [dim])
%
% Set all vectors in X along the specified dimension to be in the range
% [0,1].
%
%IN:
%   X - Array containing vectors to rescale to range [0, 1].
%   dim - Dimension along which to zero-mean X. Default: first
%         non-singleton dimension of X.
%
%OUT:
%   Y - Rescaled X.

function X = range01(X, dim)
if nargin < 2
    % Find the first non-singleton dimension
    dim = find(size(X) > 1, 1, 'first');
    if isempty(dim)
        return;
    end
end
X = bsxfun(@minus, X, min(X, [], dim));
X = bsxfun(@times, X, 1./max(X, [], dim));
end
