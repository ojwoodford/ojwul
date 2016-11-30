%IMNONMAXSUP Non-maxima suppression in an image
%
% M = imnonmaxsup(I, [radius])
%
%IN:
%   I - HxW image
%   radius - scalar value of radius (in pixels) within which to suppress
%            maxima. Default: 1.5.
%
%OUT:
%   M - HxW mask indicating which pixels are local maxima

function M = imnonmaxsup(I, radius)

% Set the default radius
if nargin < 2
    radius = 1.5;
end

% Use available non-maxima suppression function
conn = 4 + 4 * ((radius * radius) >= 2);
M = imregionalmax(I, conn);    

if radius >= 2
    % Extended supression required
    % Compute the neighbourhood
    [y, x] = ndgrid(-floor(radius):floor(radius));
    nhood = (y .* y + x .* x) <= (radius * radius);
    
    % Dilate the image
    J = imdilate(I, nhood);
    
    % Return mask of maximal points
    M = J == I & M;
end