%PROJ Project points onto the image plane
%
%   [Y, z] = proj(X)
%
% Project points into a lower dimension by dividing thorugh by the last
% dimension, e.g. pinhole camera projection.
%
%IN:
%   X - MxN input array.
%
%OUT:
%   Y - (M-1)xN output array.
%   z - 1xN array of normalizing values, z = 1./X(end,:).

function [X, z] = proj(X)
sz = size(X);
sz(1) = sz(1) - 1;
z = 1 ./ X(end,:);
X = reshape(X(1:end-1,:) .* z, sz);
sz(1) = 1;
z = reshape(z, sz);
end
