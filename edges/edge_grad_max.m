%EDGE_GRAD_MAX  Mask of gradient maxima
%
%   M = edge_grad_max(G, [thresh, [M]])
%
% Given a gradient image, computes the magnitude of locally maximal (in the
% direction of maximum gradient) gradients.
%
% IN:
%   G - HxWx2 array of gradient images in x and y directions, along third
%       dimension.
%   thresh - Scalar value determining the threshold of the second
%            derivative of gradient above which to ignore edges. Default:
%            0.
%   M - HxW logical image of pixels to look for maxima at.
%
% OUT:
%   Mo - HxWx2 gradient and gradient derivative magnitudes image, but 0 if
%        not locally maximal.
%   X - 2x(sum(M(:))) array of subpixel edgelet dimensions.

function [M, X] = edge_grad_max(G, thresh, M)

% Default values
if nargin < 3
    M = true(size(G, 1), size(G, 2));
    if nargin < 2
        thresh = 0;
    end
end

% Sample in direction of maximum gradient
GM = reshape(G, [], 2);
GM = GM(M,:);
off = bsxfun(@times, GM, 1./(max(abs(GM), [], 2) + 1e-100));
[Y, X] = find(M);
a = reshape(ojw_interp2(G, bsxfun(@plus, X, [-off(:,1) off(:,1)]), bsxfun(@plus, Y, [-off(:,2) off(:,2)]), 'l', cast(0, class(G))), [], 2, 2);

% Compute normalized edge direction
mag = normd(GM, 2);
nm = bsxfun(@times, GM, 1 ./ (mag + 1e-100));

% Compute derivatives of gradient magnitude in direction of maximum
% gradient
a = abs(sum(bsxfun(@times, a, reshape(nm, [], 1, 2)), 3));
deriv_2nd = 2 * mag - sum(a, 2);

% Output the edge gradient magnitude
M_ = M;
M = zeros(numel(M), 2);
M(M_,1) = mag;
M(M_,2) = mag;
M(M_,3) = deriv_2nd;

% Find gradient maxima above the threshold
M_(M_) = all(bsxfun(@gt, mag, a), 2) & mag > thresh & deriv_2nd > thresh;
% Mask boundary pixels
M_([1 end],:) = false;
M_(:,[1 end]) = false;

% Mask output magnitudes
M(~M_,1) = 0;
M = reshape(M, [size(G, 1), size(G, 2) 3]);

if nargout > 1
    % Sub-pixel refine the edgelets in the direction of maximum gradient
    off = 0.5 * normd(off, 2) .* diff(a, 1, 2) ./ deriv_2nd;
    X = [(X + nm(:,1) .* off)'; (Y + nm(:,2) .* off)'];
end
end
