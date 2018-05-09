%CORNERS Compute corner detector score
%
%   score = corners(I, [sigma], method)
%
%IN:
%   I - HxWxC image
%   sigma - Scalar value determing the scale of the gradient filter.
%           Default: 1.
%   method - String determing the method to use: 'harris', 'noble'
%           or 'shi-tomasi'. Default: 'shi-tomasi'.
%
%OUT:
%   score - HxW corner detector score.

function score = corners(I, varargin)
% Set the default values
method = 'shi-tomasi';
sigma = 1;
for v = varargin(:)'
    if ischar(v{1})
        method = v{1};
    elseif isscalar(v{1})
        sigma = v{1};
    else
        error('Input argument not recognized');
    end
end

% Compute the image gradient
[Ix, Iy] = imgrad(I, 'Sobel');

% Compute 1D Gaussian filter from sigma
g = max(floor((5 / 2) * sigma), 1);
g = -g:g;
g = gauss_mask(sigma, 0, g);

% Region integration
Ix2 = imfiltsep(Ix .* Ix, g', g);
Iy2 = imfiltsep(Iy .* Iy, g', g);
Ixy = imfiltsep(Ix .* Iy, g', g);

% Compute trace and determinant of structure tensor
T = Ix2 + Iy2;
D = Ix2 .* Iy2 - Ixy .* Ixy;

% Compute the corner score
switch lower(method)
    case 'harris'
        score = D - 0.04 * T .* T;
    case 'noble'
        score = D ./ (T + 1e-300);
    case 'shi-tomasi'
        score = T - sqrt(T .* T - 2 * D);
end
end 
