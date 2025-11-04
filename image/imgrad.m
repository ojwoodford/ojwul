%IMGRAD  Calculate gradients of an image in x and y directions
%
%   [Ix, Iy] = imgrad(I, sigma, [mcm])
%   G = imgrad(I, sigma, [mcm])
%
% Computes x and y gradients of an image. The implementation uses
% separable, steerable filters, hence is fast.
%
% For multi-channel images, the gradient is computed in the colour
% direction of greatest change, by defualt, using the method described in:
%    DiZenzo, S. [1986]. "A note on the gradient of a multi-image",
%    CVGIP, Vol. 33, pp. 116-125.
% Other methods can be selected.
%
% IN:
%   I - HxWxC input image.
%   filter - String denoting named filter ['Bilinear', 'Central',
%            'Prewitt', 'Sobel', 'Simoncelli'] or scalar value determing
%            the scale of the gradient filter (essentially the amount of
%            pre-smoothing to apply). If 'Simoncelli', the 7-tap values
%            from Farid, H. and Simoncelli, E. "Differentiation of Discrete
%            Multi-Dimensional Signals". IEEE Trans. Image Processing.
%            13(4): 496-508 (2004) are used.
%   mcm - multi-channel method: 'dizenzo' (default), 'norm', 'pca',
%                               'rgb2gray', 'none'.
%
% OUT:
%   Ix - HxW array of gradient images in the x direction.
%   Iy - HxW array of gradient images in the y direction.
%   G  - HxWx2 array of gradient images in x and y directions, i.e. cat(3,
%        Ix, Iy).

function [Ix, Iy] = imgrad(I, filter, mcm)

% Set the default arguments
if nargin < 3
    mcm = 'dizenzo';
end

I = double(I);
[h, w, c] = size(I);

% Convert multi-channel images to single channel, if required
if c > 1
    switch lower(mcm)
        case 'rgb2gray'
            % Standard RGB to grey conversion
            assert(c == 3, 'Image must have 3 channels');
            I = reshape(reshape(I, h*w, c) * [0.299; 0.587; 0.114], h, w);
            c = 1;
        case 'pca'
            % Compress along the first principle component of colour
            I = reshape(I, h*w, c)';
            T = pca(I);
            I = reshape(T(1,1:end-1) * I, h, w);
            c = 1;
    end
end

if ischar(filter)
    switch lower(filter)
        case 'bilinear'
            % Gradient resulting from bilinear interpolation
            g = [1 1];
            gp = [1 -1];
        case 'central'
            % Central differences
            g = 1;
            gp = [1 0 -1];
        case 'prewitt'
            g = [1 1 1];
            gp = [1 0 -1];
        case 'sobel'
            g = [1 2 1];
            gp = [1 0 -1];
        case 'bickley'
            g = [1 4 1];
            gp = [1 0 -1];
        case 'simoncelli'
            % Use 7-tap filter from:
            % Farid, H. and Simoncelli, E. "Differentiation of Discrete Multi-Dimensional Signals"
            % IEEE Trans. Image Processing. 13(4): 496-508 (2004)
            g  = [0.004711  0.069321  0.245410  0.361117  0.245410  0.069321  0.004711];
            gp = [0.018708  0.125376  0.193091  0.000000 -0.193091 -0.125376 -0.018708];
        otherwise
            error('Filter %s unrecognized', filter);
    end
    X = 1:numel(gp);
    X = X - mean(X);
else
    % Determine necessary filter support (for Gaussian).
    X = max(floor((5 / 2) * filter), 1);
    X = -X:X;
    
    % Evaluate 1D Gaussian filter (and its derivative).
    g = gauss_mask(filter, 0, X);
    gp = gauss_mask(filter, 1, X);
end

% Normalize the filters
g = g / sum(g);
gp = gp / sum(abs(gp) .* abs(X));

% Calculate image gradients
Ix = imfiltsep(I, g, gp);
Iy = imfiltsep(I, gp, g);

% Compress multi-channelled images
if c > 1
    switch lower(mcm)
        case 'dizenzo'
            % Compute maximum gradient using method in:
            % DiZenzo, S. [1986]. "A note on the gradient of a multi-image",
            % CVGIP, Vol. 33, pp. 116-125.
            
            % Compute inner product of Jacobian with itself
            J = [shiftdim(sum(Ix .* Ix, 3), -1); shiftdim(sum(Iy .* Iy, 3), -1); shiftdim(sum(Ix .* Iy, 3), -1)];
            
            % Compute the first eigen vector
            lambda = J(1,:) - J(2,:);
            lambda = 0.5 * (J(1,:) + J(2,:) + sqrt(lambda .* lambda + 4 * (J(3,:) .* J(3,:))));
            J_ = [J(2,:)-lambda; -J(3,:)];
            J = [-J(3,:); J(1,:)-lambda];
            M = sum(J .* J, 1) < sum(J_ .* J_, 1);
            J(:,M) = J_(:,M);
            J = bsxfun(@times, J, sqrt(lambda) ./ (normd(J, 1) + 1e-100));
            Ix = reshape(J(1,:), h, w);
            Iy = reshape(J(2,:), h, w);
        case 'norm'
            % Compute the gradient norm
            Ix = sign(sum(Ix, 3)) .* normd(Ix, 3);
            Iy = sign(sum(Iy, 3)) .* normd(Iy, 3);
        case 'none'
            % Leave as is
        otherwise
            error('Multi-channel method %s not recognized', mcm);
    end
end

% Concatenate if necessary
if nargout == 1
    Ix = cat(3, Ix, Iy);
end
return