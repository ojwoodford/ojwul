%EDGE_DEMO  Compute and visualize the gradient and edgels of an image
%
%   edge_demo
%   edge_demo(A, scale, thresh)
%
%IN:
%   A - HxWxC image. Default: Use peppers.png
%   scale - sigma of Gaussian blur to apply during edge detection.
%           Default: 0.7.
%   thresh - Threshold on edgel suppression. If negative, magnitude is
%            interpreted as the proportion of pixels to return as edgels.
%            Default: -0.05, i.e. keep top 5% of pixels, by edge strength.

function edge_demo(A, scale, thresh)
% Set default arguments
if nargin < 3
    thresh = -0.05; % Negative numbers indicate ratio of edges to keep
    if nargin < 2
        scale = 0.7;
        if nargin < 1
            A = imread('peppers.png');
        end
    end
end

% Compute the gradient image
tic;
G = imgrad(A, scale, 'dizenzo');
t1 = toc;

% Compute a reasonable edgel threshold
if thresh < 0
    s = sort(reshape(sum(G .* G, 3), [], 1), 'descend');
    thresh = sqrt(s(ceil(end*-thresh)));
    fprintf('Edgel threshold: %g\n', thresh);
end

% Compute the edgels
tic;
X = compute_edgels(G, thresh);
t2 = toc;

% Visualization
fprintf('Computing gradients: %gs\nExtracting edgels: %gs\n', t1, t2);
clf;
% Render the gradient image
sc(G, 'flow', [0 thresh*2]);
% Render the edgels
hold on
c = 0.5 * cos(reshape(X(3,:), 1, 1, []));
s = 0.5 * sin(reshape(X(3,:), 1, 1, []));
X = bsxfun(@plus, reshape(X(1:2,:), 2, 1, []), [-s s; c -c]);
plot(squeeze(X(1,:,:)), squeeze(X(2,:,:)), 'r-', 'LineWidth', 2);
hold off