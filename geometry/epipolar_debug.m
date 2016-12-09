%EPIPOLAR_DEBUG  Show epipolar geometry between two images
%
%   handles = epipolar_debug(F, X, handles)
%
% Given a pair of correspondences and their fundamental matrix, or two
% projection matrices from world to image coordinates, this function
% renders the epipolar geometry for each correspondence in the two images.
%
%IN:
%   F - 3x3 fundamental matrix, 3x4 projection matrix or 3x4x2 set of
%       projection matrices.
%   X - 4xN set of 2D image correspondences, in calibrated coordinates.
%   handles - Set of graphics handles generated in previous call to
%             epipolar_debug.
%
%OUT:
%   handles - Set of graphics handles for updates.

function handles = epipolar_debug(F, X, handles)
% Convert projection matrices to a fundamental matrx
if isequal(size(F), [3 4])
    F = cat(3, eye(3, 4), F);
end
if isequal(size(F), [3 4 2])
    F = F_from_P(F);
end

if nargin < 3
    % Set up the plots
    handles{1} = 2 * max(abs(reshape(X, 2, [])), [], 2);
    figure(6296);
    clf reset;
    subplot(121);
    handles = [handles; render_image(F, X, handles{1})];
    subplot(122);
    handles = [handles; render_image(F', X([3:4 1:2],:), handles{1})];
else
    % Update the plots
    h = reshape(handles(2:end), [], 2);
    render_image(F, X, handles{1}, h(:,1));
    render_image(F', X([3:4 1:2],:), handles{1}, h(:,1));
end
end

function handles = render_image(F, X, sz, handles)
% Compute the epipolar lines from the other image
n = size(X, 2);
L = F' * homg(X(3:4,:));
x = repmat([-sz(1); sz(1)], 1, n);
y = repmat([-sz(2); sz(2)], 1, n);
M = L(1,:) > L(2,:);
y(:,M) = (-L(3,M) - L(1,M) .* x(:,M)) ./ L(2,M);
M = ~M;
x(:,M) = (-L(3,M) - L(2,M) .* x(:,M)) ./ L(1,M);
% Render
if nargin < 4
    % Plot the point to line distances
    handles{1} = plot(X(1,:)', X(2,:)', 'b.', 'MarkerSize', 20);
    % Plot the epipolar lines
    handles{2} = plot(x, y, 'r-');
    % Set the frame size
    xlim([-sz(1) sz(1)]);
    ylim([-sz(2) sz(2)]);
else
    % Updates
    set(handles{1}, 'XData', X(1,:)', 'YData', X(2,:)');
    set(handles{2}, 'XData', x, 'YData', y);
end
end
   