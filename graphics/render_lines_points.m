%RENDER_LINES_POINTS Renders a set of coloured points and/or lines
%
%   h = render_line_points(X, [C, [options]])
%
% Renders a set of coloured 2D or 3D points or lines, either as a set of
% lineseries objects (quantized colours, fast rendering) or as a patch
% object (exact colours, slow rendering), depending on the number of
% colours requested.
%
%IN:
%   X - MxNxP set of P M-D points (if N==1) or lines with N-1 segments,
%       where M is 2 or 3.
%   C - 3xP set of RGB colours for each point/line, or a single colorspec.
%       Default: 'b'.
%   options - string value pairs:
%      NumColors - scalar indicating the number of colour quantization
%                  levels (0 for no quantization). Default: 0.
%      Marker - character representing the marker to use. Default: '.' for
%               points, 'none' for lines.
%      MarkerSize - scalar marker size. Default: 6.
%      LineStyle - line style string. Default: '-'.
%      LineWidth - scalar line width. Default: 0.5.
%      Alpha - scalar line transparency. Default: 1.
%
%OUT:
%   h - handle(s) to created graphics object(s).

function h = render_lines_points(X, C, varargin)
% Set defaults
if size(X, 2) == 1
    opts.Marker = '.';
else
    opts.Marker = 'none';
end
opts.MarkerSize = 6;
opts.LineStyle = '-';
opts.LineWidth = 0.5;
opts.NumColors = 0;
opts.Alpha = 1;
if nargin < 2
    C = 'b';
else
    % Parse input options
    opts = vgg_argparse(opts, varargin);
end
% Create the options parameters
alpha = opts.Alpha;
if alpha < 1
    num_colors = 0;
else
    num_colors = opts.NumColors;
end
if size(X, 2) == 1
    opts.LineStyle = 'none';
end
opts = {'Marker', opts.Marker, 'MarkerSize', opts.MarkerSize, 'LineWidth', opts.LineWidth, 'LineStyle', opts.LineStyle};

X(:,end+1,:) = NaN;
col = @(X) X(:);
if numel(C) <= 3 && alpha == 1
    % One colour, so just plot as a line series
    if size(X, 1) < 3
        h = plot(col(X(1,:,:)), col(X(2,:,:)), 'b-', 'Color', C(:)', opts{:});
    else
        h = plot3(col(X(1,:,:)), col(X(2,:,:)), col(X(3,:,:)), 'b-', 'Color', C(:)', opts{:});
    end
    return;
end
% Colour per point
C = C';
if num_colors < 1
    % Plot as a patch object: points are exact colours, but the plot is
    % slow to render and manipulate
    if isscalar(C)
        C = repmat(rem(floor((strfind('kbgcrmyw', C) - 1) * [0.25 0.5 1]), 2), [size(X, 3) 1]);
    end
    C = repmat(reshape(C, 1, [], 3), [size(X, 2) 1 1]);
    if size(X, 1) < 3
        h = patch(shiftdim(X(1,:,:), 1), shiftdim(X(2,:,:), 1), C, 'FaceColor', 'none', 'EdgeColor', 'flat', 'MarkerEdgeColor', 'flat', opts{:}, 'EdgeAlpha', alpha);
    else
        h = patch(shiftdim(X(1,:,:), 1), shiftdim(X(2,:,:), 1), shiftdim(X(3,:,:), 1), C, 'FaceColor', 'none', 'EdgeColor', 'flat', 'MarkerEdgeColor', 'flat', opts{:}, 'EdgeAlpha', alpha);
    end
else
    % Plot as a set of single colour line objects: colours are quantized
    % but the plot is faster to update
    [I, C] = rgb2ind(reshape(C, 1, [], 3), num_colors, 'nodither');
    num_colors = size(C, 1);
    I = accumarray(I(:)+1, (1:size(X, 3))', [num_colors 1], @(I) {I});
    tf = ishold();
    hold on
    for a = num_colors:-1:1
        if ~isempty(I{a})
            if size(X, 1) < 3
                h(a) = plot(col(X(1,:,I{a})), col(X(2,:,I{a})), 'b-', 'Color', C(a,:), opts{:});
            else
                h(a) = plot3(col(X(1,:,I{a})), col(X(2,:,I{a})), col(X(3,:,I{a})), 'b-', 'Color', C(a,:), opts{:});
            end
        end
    end
    if ~tf
        hold off
    end
end
