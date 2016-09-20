%RENDER_CONFUSION Render a confusion matrix
%
%   render_confusion(M[, class_names])
%
% Renders a confusion matrix into the current axes.
%
%IN:
%   M - NxN confusion matrix.
%   class_names - Nx1 cell array of class name strings. Default: {} (not
%                 shown).

function render_confusion(M, class_names, varargin)

% Default optional arguments
opts.cmap = '-gray'; 
%opts.cmap = '-hot*'; % Alternative option
opts.totals = false;
opts.labels = true;
opts.font_size = 0.4;  % Proportion of box height
opts.font = 'Times';
opts = vgg_argparse(opts, varargin);

% Create a totals row at the bottom
if opts.totals
    M(end+1,:) = sum(M, 1);
end

% Create the background image
im = sc(kron(bsxfun(@times, M, 1./sum(M, 2)), ones(8)), opts.cmap, [0 1]);
image(im);
im = sc(im(1:8:end,1:8:end,:), 'rgb2gray');
im = im(:,:,1);

% Write the values
[h, w] = size(M);
color = {'k', 'w'};
hold on
for x = 1:w
    for y = 1:h
        text(x*8-3.5, y*8-3.5, sprintf('%.2g', M(y,x)), 'FontName', opts.font, 'Units', 'Data', 'FontUnits', 'Normalized', 'FontSize', opts.font_size/h, 'Color', color{1+(im(y,x)<0.5)}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
hold off

% Write labels
if opts.labels
    set(xlabel('Output class'), 'FontName', opts.font, 'FontUnits', 'Normalized', 'FontSize', 1.2*opts.font_size/h);
    set(ylabel('Ground truth class'), 'FontName', opts.font, 'FontUnits', 'Normalized', 'FontSize', 1.2*opts.font_size/h);
end
if opts.totals
    text(w*8+2, h*8-3.5, 'TOTAL', 'FontName', opts.font, 'Units', 'Data', 'FontUnits', 'Normalized', 'FontSize', opts.font_size/h, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

% Sort out the axes
axis equal tight
set(gca, 'XTick', [], 'YTick', []);

if nargin < 2 || isempty(class_names)
    return
end

% Write the class names
hold on
for x = 1:numel(class_names)
    if opts.labels
        % Right edge
        text(w*8+2, x*8-3.5, class_names{x}, 'FontName', opts.font, 'Rotation', 45, 'Units', 'Data', 'FontUnits', 'Normalized', 'FontSize', opts.font_size/h, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        % Top edge
        text(x*8-3.5, -2, class_names{x}, 'FontName', opts.font, 'Rotation', 45, 'Units', 'Data', 'FontUnits', 'Normalized', 'FontSize', opts.font_size/h, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    else
        % Left edge
        text(-2, x*8-3.5, class_names{x}, 'FontName', opts.font, 'Rotation', 45, 'Units', 'Data', 'FontUnits', 'Normalized', 'FontSize', opts.font_size/h, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        % Bottom edge
        text(x*8-3.5, h*8+2, class_names{x}, 'FontName', opts.font, 'Rotation', 45, 'Units', 'Data', 'FontUnits', 'Normalized', 'FontSize', opts.font_size/h, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
end
hold off
