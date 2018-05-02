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
radius2 = radius * radius;

% Find 4-connected maxima (must be greater than neighbours)
M = I(2:end-1,2:end-1);
M = M > I(1:end-2,2:end-1) & ...
    M > I(3:end,2:end-1)   & ...
    M > I(2:end-1,1:end-2) & ...
    M > I(2:end-1,3:end);
M = [false(1, size(M, 2)+2); false(size(M, 1), 1) M false(size(M, 1), 1); false(1, size(M, 2)+2)];

if radius2 >= 2
    % Extended supression required
    % Compute the neighbourhood
    [y, x] = ndgrid(-floor(radius):floor(radius));
    nhood = (y .* y + x .* x) <= radius2;
    
    % Dilate the image
    try
        J = imdilate(I, nhood);
    catch
        x = x(nhood);
        y = y(nhood);
        [Y, X] = ndgrid(1:size(I, 1), 1:size(I, 2));
        J = reshape(max(ojw_interp2(I, bsxfun(@plus, X(:)', x(:)), bsxfun(@plus, Y(:)', y(:)), 'n', cast(0, class(I)))), size(I)); 
    end
    
    % Return mask of maximal points
    M = J == I & M;
end
end