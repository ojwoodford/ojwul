%EPIPOLAR_DEMO  Interactively show epipolar geometry between two images
%
%   epipolar_demo(ims, F)
%
% Given a pair of images and their fundamental matrix, or two projection
% matrices from world to image coordinates, this function renders the
% epipolar geometry of a given pixel from one image in the second image.
% Lines can be stored by clicking on pixels.
%
% IN:
%   ims - 2x1 cell array of input images
%   F - 3x3 fundamental matrix or 3x4x2 set of projection matrices.

function epipolar_demo(ims, F)

% Convert projection matrices to a fundamental matrx
if isequal(size(F), [3 4 2])
    F = F_from_P(F);
end
state.F = {F, F'};

% Initialize the plot
clf;
nn = [NaN; NaN];
for a = 2:-1:1
    % Render the image
    state.im(a).hAx = axes('Position', [0.01+(a-1)*0.495 0.01 0.485 0.98]);
    state.im(a).sz = size(ims{a});
    image(ims{a});
    colormap(gray(256));
    axis equal off
    xlim([0.5 state.im(a).sz(2)+0.5]);
    ylim([0.5 state.im(a).sz(1)+0.5]);
    % Render the line
    hold on
    state.im(a).hL = plot(nn, nn, 'r-');
    hold off
end

set(gcf(), 'Interruptible', 'off', 'BusyAction', 'cancel', 'WindowButtonUpFcn', @mouseup_callback, 'WindowButtonMotionFcn', @mouse_callback, 'UserData', state);
end

function mouse_callback(fig, evnt)
state = get(fig, 'UserData');
nn = [NaN; NaN];
for a = 1:2
    % Get the mouse location
    y = get(state.im(a).hAx, 'CurrentPoint');
    x = round(y(1,1));
    y = round(y(1,2));
    % Check if it is over the image
    if x < 1 || x > state.im(a).sz(2) || y < 1 || y > state.im(a).sz(1)
        % Clear the line in the other image
        X = nn;
        Y = nn;
    else
        % Plot the epipolar line in the other image
        [X, Y] = compute_line(state, a, x, y);
    end
    set(state.im(3-a).hL, 'XData', X, 'YData', Y);
end
drawnow();
end

function mouseup_callback(fig, evnt)
state = get(fig, 'UserData');
% Render the point and line in a random colour
col = rand(1, 3);
for a = 1:2
    % Get the mouse location
    y = get(state.im(a).hAx, 'CurrentPoint');
    x = round(y(1,1));
    y = round(y(1,2));
    % Check if it is over the image
    if x < 1 || x > state.im(a).sz(2) || y < 1 || y > state.im(a).sz(1)
        continue;
    end
    % Plot the epipolar line in the other image
    [X, Y] = compute_line(state, a, x, y);
    set(fig, 'CurrentAxes', state.im(a).hAx);
    hold on
    plot(x, y, 'r+', 'Color', col);
    hold off
    set(fig, 'CurrentAxes', state.im(3-a).hAx);
    hold on
    plot(X, Y, 'r-', 'Color', col);
    hold off
end
drawnow;
end

function [X, Y] = compute_line(state, a, x, y)
L = state.F{a} * [x; y; 1];
if L(1) > L(2)
    X = [0.5; state.im(3-a).sz(2)+0.5];
    Y = (-L(3) - L(1) * X) / L(2);
else
    Y = [0.5; state.im(3-a).sz(1)+0.5];
    X = (-L(3) - L(2) * Y) / L(1);
end
end

   