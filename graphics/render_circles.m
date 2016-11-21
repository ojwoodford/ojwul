%RENDER_CIRCLES Draw a set of circles
%
%   h = render_circles(X, nsides, col)
%
%IN:
%   X - 3xN array of circles defined by [center_x; center_y; radius].
%   nsides - Scalar indicating the number of sides to use for each polygon.
%            Default: 64.
%   col - Colour to use for the circles. Default: 'b'.

function h = render_circles(X, nsides, color)
    
% Default parameters
if nargin < 2
    nsides = 64;
end
nsides = max(round(nsides), 3); % Make sure it is an integer >= 3
if nargin < 3
    color = [0 0 1];
end

% Compute the coordinates
X = reshape(X, 3, []);
radii = X(3,:);
poly_sides = linspace(0, 2 * pi, nsides)';
Y = bsxfun(@plus, bsxfun(@times, sin(poly_sides), radii), X(2,:));
X = bsxfun(@plus, bsxfun(@times, cos(poly_sides), radii), X(1,:));

% Expand the color
if size(color, 1) ~= size(X, 2)
    color = repmat(color, [size(X, 2) 1]);
end

% Render using patch (so can have different coloured circles)
h = patch(X, Y, repmat(shiftdim(color, -1), [size(X, 1) 1 1]), 'FaceColor', 'none', 'EdgeColor', 'flat');
end
