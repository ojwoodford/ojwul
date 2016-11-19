%IMNORM Spatially local image normalization
%
%   B = imnorm(A, sigma)
%   B = imnorm(A, [szy szx])
%   B = imnorm(..., nonlin)
%
% Apply a local normalization operator to an image, either with a Gaussian
% or window average weighting.
%
%IN:
%   A - HxWxC input image.
%   sigma - scalar indicating the standard deviation of the Gaussian
%           weighting to apply, in pixels.
%   [szy szx] - 1x2 window size for window average weighting.
%   nonlin - a handle to an elementwise non-linear function to be applied
%            to the denominator prior to normalization.
%
%OUT:
%   B - HxWxC normalized image.

function A = imnorm(A, sigma, nonlin)
% Compute the separable filters
if isscalar(sigma)
    fx = max(floor(5/2 * sigma), 1);
    fx = gauss_mask(sigma, 0, -fx:fx);
    fx = fx / sum(fx);
    fy = fx;
elseif isequal(size(sigma), [1 2]) && isequal(sigma, round(sigma))
    fy = zeros(1, sigma(1)) + 1 / sigma(1);
    fx = zeros(1, sigma(2)) + 1 / sigma(2);
else
    error('sigma not recognized');
end

% Filter the image
A = double(A);
B = imfiltsep(A, fy, fx);

% Apply the non-linear function
if nargin > 2
    B = nonlin(B);
end

% Normalize the image
A = A ./ B;
end
