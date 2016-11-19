%EDGE_SCALE  Compute an edge scale image
%
%   [S, G, E] = edge_scale(I, [max_octaves, [levels_per_octave, [thresh]]])
%
% Output the edge scale of each pixel, and a list of edgels, if requested.
% The scale of an edge is determined by the amount of image smoothing
% required before it stops being a locally maximal gradient.
%
% IN:
%   I - HxWxC input image.
%   max_octaves - integer indicating the number of scale octaves to search
%                 over, from 0 (original image only) to Inf (coarsest scale
%                 possible). Default: Inf.
%   levels_per_octave - positive integer indicating the number scale levels
%                       to have per octave (halving of image size).
%                       Default: 1.
%
% OUT:
%   S - HxW edge scale image, indicating for each pixel the coarsest octave
%       at which that pixel is an edge (starting from 0), or -1 if not an
%       edge.
%   G - HxWx2 gradient image at the finest scale.
%   E - 5xN list of edgels, each column giving [x y angle mag scale] of the
%       edgel centre.

function [S, G, E] = edge_scale(I, max_octaves, levels_per_octave, thresh)

% Default values
if nargin < 4
    thresh = 0;
    if nargin < 3
        levels_per_octave = 1;
        if nargin < 2
            max_octaves = Inf;
        end
    end
end
[h, w, c] = size(I);
max_octaves = min(max_octaves, floor(log2(min(h, w))) - 3);
max_octaves = max_octaves * levels_per_octave - 1;

% Prefilter the image
sigma = 1.0;
if sigma > 0.5
    I = filter(I, sqrt(sigma * sigma - 0.25));
end

% Compute the gradient image at the coarsest scale
grad = @(I) imgrad(I, 0, 'norm');
[Ix, Iy] = grad(I);
G = cat(3, Ix, Iy);

if nargout < 3
    % Compute the gradient maxima and subpixel locations
    M = edge_grad_max(G, thresh);
    M = M(:,:,1) ~= 0;
else
    % Compute the gradient maxima and subpixel locations
    [M, X] = edge_grad_max(G, thresh);
    mag = M(:,:,1);
    M = mag ~= 0;
    mag = mag(M);
    
    % Compute the edge angles, in radians
    angle = atan2(Iy(M), Ix(M));
    
    % Output only those edgelets which meet the criteria
    E = [X(:,M); angle'; mag'];
    Morig = M;
end
clear Ix Iy

% Initialize the scale image
S = (M - 1) * levels_per_octave;

% Now go over each further octave
scale = 0;
octave = 0;
Mlow = M;
while any(Mlow(:))
    % Increase the scale
    scale = scale + 1;
    if scale > max_octaves
        break;
    end
    
    % Compute the filter
    sigma_new = sigma * (2 ^ (1 / levels_per_octave));
    g = gauss_mask(sqrt(sigma_new * sigma_new - sigma * sigma));
    sigma = sigma_new;
    
    % Filter the image
    I = filter(I, g);
    
    % Subsample if necessary
    if mod(scale, levels_per_octave) == 0
        I = I(1:2:end,1:2:end,:);
        Mlow = conv2([1 1 1], [1; 1; 1], single(Mlow), 'same') > 0;
        Mlow = Mlow(1:2:end,1:2:end); 
        sigma = sigma / 2;
        octave = octave + 1;
    end
    
    % Find the gradient maxima
    M_ = edge_grad_max(grad(I), 0, Mlow);
    M_ = M_(:,:,1) ~= 0;
    Mlow = Mlow & M_;
    
    % Upscale to original size
    M_(end+1,end+1) = false;
    M_ = single(M_);
    for a = 1:octave;
        M_ = interp2(M_);
    end
    
    % Update the mask
    M_ = M_(1:h,1:w) ~= 0;
    M = M & M_;
    
    % Update the scales
    S(M) = scale;
end

% Rescale the scale
S = S / levels_per_octave;

% Add scale onto the edgelet descriptions
if nargout > 1
    E = [E; S(Morig)'];
end
end

function I = filter(I, g)
g = g / sum(g);
I = imfiltsep(I, g, g);
end