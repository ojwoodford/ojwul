%CONVERT2GRAY Convert an RGB image to grayscale
%
%   B = convert2gray(A)
%
%IN:
%   A - HxWxC input image, where C = 3 (RGB) or 1 (already grayscale).
%
%OUT:
%   B - HxW grayscale output image, of the same class as A.

function im = convert2gray(im)
if size(im, 3) == 3
    im = cast(reshape(double(reshape(im, [], 3)) * [0.299; 0.587; 0.114], size(im, 1), size(im, 2)), 'like', im);
end
end
