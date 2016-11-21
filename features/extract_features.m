%EXTRACT_FEATURES Extract interest points from a detector score image
%
%   [X, s] = extract_features(score, [radius, [thresh, [subpixel]]])
%
% Finds the local maxima in the input detector score image.
%
%IN:
%   score - HxW interest point detector score.
%   radius - Scalar indicating the radius of non-maxima suppression to
%            apply. Default: 1.5.
%   thresh - Threshold on corner suppression. If negative, magnitude is
%            interpreted as the proportion of strongest corners to keep.
%            Default: -0.5, i.e. keep top 50% of corners.
%   subpixel - Boolean indicating whether subpixel refinement
%   
%OUT:
%   X - 2xN matrix of interest point locations.
%   s - 1xN vector of interest point scores.

function [X, score] = extract_features(score, radius, thresh, subpixel)
% Set the default values
if nargin < 4
    subpixel = true;
    if nargin < 3
        radius = 1.5;
        if nargin < 2
            thresh = -0.5;
        end
    end
end

% Find the local maxima
M = imnonmaxsup(score, radius);

% Reject boundary maxima
M([1 end],:) = false;
M(:,[1 end]) = false;

% Keep only those above the threshold
if thresh < 0
    s = sort(score(M), 'descend');
    thresh = s(ceil(end*-thresh));
end
M = M & score >= thresh;

% Compute the indices
[y, x] = find(M);
X = [y(:)'; x(:)'];

if subpixel
    % Sub-pixel refinement
    X = X + refine_subpixel(score, M);
end

% Put x coordinate first
X = X([2 1],:);

% Output the scores
if nargout > 1
    score = score(M);
end