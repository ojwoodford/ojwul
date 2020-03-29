%IMDILATE_ Dilate image
%
%   B = imdilate_(A, nhood)
%
% This function duplicates behaviour of the IMDILATE function (without
% needing the Image Processing Toolbox) for a subset of inputs. If IMDILATE
% is available, it is used.
%
%IN:
%   A - HxW logical or grayscale image.
%   nhood - PxQ matrix of 1s and 0s indicating the neighbourhood of a pixel
%           to be dilated.
%
%OUT:
%   B - HxW dilated image.

function B = imdilate_(A, nhood)

% Use the IPT if available
if license('test', 'image_toolbox')
    B = imdilate(A, nhood);
    return;
end

% Cache indices for the dilation
persistent cache
if isempty(cache) || ~isequal(nhood, cache.nhood) || ~isequal(size(A), cache.sz)
    cache.nhood = nhood;
    cache.sz = size(A);
    [y, x] = find(nhood);
    y = int32(y - ceil(size(nhood, 1) * 0.5 + 0.5));
    x = int32(x - ceil(size(nhood, 2) * 0.5 + 0.5)) - int32(1);
    [Y, X] = ndgrid(int32(1:cache.sz(1)), int32(1:cache.sz(2)));
    X = min(max(X(:)' + x(:), int32(0)), cache.sz(2)-1);
    Y = min(max(Y(:)' + y(:), int32(1)), cache.sz(1));
    cache.ind = X * int32(cache.sz(1)) + Y;
end

% Perform the dilation
B = reshape(max(A(cache.ind), [], 1), cache.sz);
end