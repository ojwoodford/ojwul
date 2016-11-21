%HARRIS Harris corner detector
%
%   score = harris(I, [sigma])
%
%IN:
%   I - HxWxC image
%   sigma - Scalar value determing the scale of the gradient filter.
%           Default: 1.
%
%OUT:
%   score - HxW corner detector score.

function score = harris(I, sigma)
% Set the default values
if nargin < 2
    sigma = 1;
end

% Compute the image gradient
[Ix, Iy] = imgrad(I, 0);

% Compute 1D Gaussian filter from sigma
g = max(floor((5 / 2) * sigma), 1);
g = -g:g;
g = gauss_mask(sigma, 0, g);

% Region integration
Ix2 = imfiltsep(Ix .* Ix, g', g);
Iy2 = imfiltsep(Iy .* Iy, g', g);
Ixy = imfiltsep(Ix .* Iy, g', g);

% Compute the corner score (Peter Kovesi's)
score = (Ix2 .* Iy2 - Ixy .* Ixy) ./ (Ix2 + Iy2 + 1e-300);
end 
