%IMWARP  Warp an image according to a homography
%
%   im = imwarp(im, H)
%
%IN:
%   im - HxWxC image
%   H - 3x3 homography matrix from source to target
%
%OUT:
%   im - HxWxC resampled output image

function im = imwarp(im, H)
% Compute the coordinates to sample at
X = proj(H \ homg(flipud(ndgrid_cols(1:size(im, 1), 1:size(im, 2)))));
% Resample the image
im = reshape(ojw_interp2(im, X(1,:), X(2,:), '6', 0), size(im));
end
