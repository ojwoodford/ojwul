function h = render_frame(A, P, scale, color)
%RENDER_FRAME Renders an image, A, and a camera in the position given by P

scale = diag([abs(scale) abs(scale) scale]);

% Render the camera
if nargin < 4
    color = [0.5 0.5 0.5];
end
[K R t] = KR_from_P(P);
h = render_camera([0.1 * scale * R', -t], color);

if ~isempty(A)
    % Render the frame
    X = (size(A, 2) / K(1)) * [-0.5 0.5];
    Y = (size(A, 1) / K(5)) * [-0.5 0.5];
    [Y X] = ndgrid(Y, X);
    W = reshape(([scale * R', -t] * [X(:)'; Y(:)'; ones(2, numel(X))])', 2, 2, 3);
    h(2) = surface(W(:,:,1), W(:,:,2), W(:,:,3), im2double(A), 'facecolor', 'texturemap', 'edgecolor', 'none');
end