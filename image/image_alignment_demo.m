%IMAGE_ALIGNMENT_DEMO Demonstrates image alignment
%
%   image_alignment_demo()
%
% This function demonstrates image alignment via a homography, using
% various techniques:
%   - Lie algebra: a Lie representation of homography updates.
%   - Auto differentiation: automatically compute the Jacobian of residuals.
%   - Gauss-Newton: gradient descent with automatic step size computation.
%   - Analytic image gradients: fast image & gradient sampling using ojw_interp2.

function image_alignment_demo(A)
if nargin < 1
    A = imread('peppers.png');
end
sz = size(A);
render();

% Create a target image
[y, x] = ndgrid(linspace(sz(1)/2-50, sz(1)/2+50, 100), linspace(sz(2)/2-50, sz(2)/2+50, 100));
ref = ojw_interp2(A, x, y); % Did you call "mex ojw_interp2.cpp" first?

% Create a small homography
Hgt = eye(3) + randn(3) * 0.0002;
X = Hgt * [x(:)'; y(:)'; ones(1, numel(x))];

% Use SL(3) representation of a homography
sl3 = lie('sl3');

% Create a cost function
H = eye(3);
for a = 1:100
    % Create a delta vector of 8 variables
    dH = autodiff(zeros(8, 1)); % The magic line!!!!
    
    % Compute the homgraphy update
    dH = exp(sl3, dH);
    
    % Apply the homography
    Y = (H * dH) * X;
    Y = bsxfun(@times, Y(1:2,:), Y(3,:));
    
    % Sample the image
    tgt = ojw_interp2(A, Y(1,:), Y(2,:));
    tgt = reshape(tgt, size(ref));

    % Compute the residuals
    r = tgt - ref;
    
    % Read out the *automagically* computed gradient of the residuals
    J = reshape(grad(r), 8, []);
    r = double(r);
    
    % Compute the gauss-newton step
    step = J' \ r(:);
    
    % Visualization
    render(ref, r, tgt);
    
    % Check for convergence
    if norm(step) < 1e-12
        break;
    end
    
    % Apply the step
    H = H * exp(sl3, -step);
end
end

function render(ref, r, tgt)
persistent handles
if nargin == 0
    handles = [];
    return;
end
cost = r(:)' * r(:);
r = log(1 + sum(r .* r, 3));
tgt = double(tgt) / 255;
try
    handles.X = [handles.X handles.X(end)+1];
    handles.Y = [handles.Y cost];
    figure(handles.fig);
    set(handles.diff_im, 'CData', r);
    set(handles.tgt_im, 'CData', tgt);
    set(handles.plot, 'XData', handles.X, 'YData', handles.Y);
catch
    handles.fig = figure(5680);
    clf;
    set(handles.fig, 'Color', 'w');
    axes('OuterPosition', [0 0.5 1/3 0.5]);
    image(ref/255);
    title 'Reference image'
    axis equal off
    axes('OuterPosition', [1/3 0.5 1/3 0.5]);
    handles.diff_im = imagesc(r);
    title 'Pixel difference magnitude'
    axis equal off
    axes('OuterPosition', [2/3 0.5 1/3 0.5]);
    handles.tgt_im = image(tgt);
    title 'Target image'
    axis equal off
    axes('OuterPosition', [0 0 1 0.5]);
    handles.X = 1;
    handles.Y = cost;
    handles.plot = plot(handles.X, handles.Y);
    xlim([1 100]);
    ylim([0 cost*1.05]);
    ylabel Cost
    xlabel Iteration
end
drawnow;
end
