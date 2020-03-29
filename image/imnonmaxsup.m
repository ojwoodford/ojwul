%IMNONMAXSUP Non-maxima suppression in an image
%
% M = imnonmaxsup(I, [radius, [mask_val]])
%
%IN:
%   I - HxW image
%   radius - scalar value of radius (in pixels) within which to suppress
%            maxima. Default: 1.5.
%   mask_val - scalar value which is used to mask previous maxima prior to
%              iterating the search for maxima. Default: output single
%              iteration maxima only.
%
%OUT:
%   M - HxW mask indicating which pixels are local maxima

function M = imnonmaxsup(I, radius, mask_val)

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

if radius2 >= 2 || nargin > 2
    % Extended supression required
    % Compute the neighbourhood
    [y, x] = ndgrid(-floor(radius):floor(radius));
    nhood = (y .* y + x .* x) <= radius2;
    
    % Dilate the image
    J = imdilate_(I, nhood);
    
    % Return mask of maximal points
    M = J == I & M;
    
    if nargin > 2
        % Iterate the search
        I(imdilate_(M, nhood)) = mask_val;
        while any(I(:))
            J = imdilate_(I, nhood);
            J = (J == I) & (I ~= mask_val);
            I(imdilate_(J, nhood)) = mask_val;
            M = M | J;
        end
    end
end
end