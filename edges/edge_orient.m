%EDGE_ORIENT  Calculate edge filter responses at any orientation
%
%   [R, O, G] = edge_orient(I, sigma, [mcm])
%
% Returns a handle to a function that generates an edge filter response
% image using angularly adaptive filtering, and also outputs the angle that
% gives maximal edge response. The implementation uses separable, steerable
% filters, hence is fast.
%
% For multi-channel images, the gradient is computed in the colour
% direction of greatest change.
%
% IN:
%   I - HxWxC input image.
%   sigma - Scalar value determing the scale of the edge filter
%           (essentially the amount of pre-smoothing to apply).
%   mcm - multi-channel method: 'dizenzo' (default), 'norm', 'pca',
%                               'rgb2gray'.
%
% OUT:
%   R - @(theta) anonymous function, which, when given a scalar angle (in
%       radians), or HxW angle array, returns an HxW edge filter response
%       image for the angle(s) given. R(O) is the maximal edge filter
%       response for each pixel.
%   O - HxW edge orientation image, returning the angle, in radians, of the
%       maximal filter response at each pixel.
%   G - HxWx2 array of gradient images in x and y directions, along third
%       dimension.

function [R, O, G] = edge_orient(I, sigma, mcm)

% Set the default method
if nargin < 3
    mcm = 'dizenzo';
end

% Compute the image gradients
[Ix, Iy] = imgrad(I, sigma, mcm);

% Output function that generates response at any angle
R = @(theta) cos(theta) .* Ix + sin(theta) .* Iy;

if nargout > 1
    % Evaluate the orientation, in radians
    O = atan2(Iy, Ix);
    if nargout > 2
        G = cat(3, Ix, Iy);
    end
end
return