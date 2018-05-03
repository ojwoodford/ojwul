%IMNORM Spatially local image normalization
%
%   B = imnorm(A, sigma, noise_variance)
%   B = imnorm(A, [szy szx], noise_variance)
%
% Apply a local normalization operator (subtracting the mean and
% normalizing the variance) to an image, either with a Gaussian or window
% average weighting.
%
%IN:
%   A - HxWxC input image.
%   sigma - scalar indicating the standard deviation of the Gaussian
%           weighting to apply, in pixels.
%   [szy szx] - 1x2 window size for window average weighting.
%   noise_variance - scalar variance of image noise.
%
%OUT:
%   B - HxWxC normalized image.

function A = imnorm(A, sigma, noise_variance)
% Compute the separable filters
if isscalar(sigma)
    fx = max(floor(5/2 * sigma), 1);
    fx = gauss_mask(sigma, 0, -fx:fx);
    fx = fx / sum(fx);
    fy = fx;
elseif isequal(size(sigma), [1 2]) && isequal(sigma, round(sigma))
    fy = repmat(1 / sigma(1), sigma(1), 1);
    fx = repmat(1 / sigma(2), 1, sigma(2));
else
    error('sigma not recognized');
end

% Filter the image
A = double(A);

% Subtract the mean
A = A - imfiltsep(A, fy, fx);

% Divide by the standard deviation
noise_variance = max(noise_variance, 1e-100);
A = A ./ sqrt(imfiltsep(A .* A, fy, fx) + noise_variance);
end
