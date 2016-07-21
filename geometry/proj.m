%PROJ Project points onto the image plane
%
%   Y = proj(X)
%
% Project points into a lower dimension by dividing thorugh by the last
% dimension, e.g. pinhole camera projection.
%
%IN:
%   X - MxN input array.
%
%OUT:
%   Y - (M-1)xN output array.

function X = proj(X)
sz = size(X);
sz(1) = sz(1) - 1;
X = reshape(bsxfun(@times, X(1:end-1,:), 1 ./ X(end,:)), sz);
end
