%RENDER_DEPTHMAP Renders a coloured depthmap as a point cloud
%
%   h = render_depthmap(P, Z, [im, [num_colors, [marker, [marker_size]]]])
%
% Renders a coloured depthmap in world coordinates.
%
%IN:
%   P - 3x4 projection matrix from world to image coordinates.
%   Z - MxN depthmap.
%   im - MxNx3 RGB image. Default: [] (points coloured blue);
%   num_colors - scalar indicating the number of colour quantization levels
%                (0 for no quantization). Default: 0.
%   marker - character representing the marker to use. Default: '.'.
%   marker_size - scalar marker size. Default: 6.
%
%OUT:
%   h - handle(s) to created graphics object(s).

function h = render_depthmap(P, Z, im, num_colors, marker, marker_size)
% Set defaults
if nargin < 6
    marker_size = 6;
    if nargin < 5
        marker = '.';
        if nargin < 4
            num_colors = 256;
            if nargin < 3
                im = [];
            end
        end
    end
end
if isempty(im)
    im = 'b';
else
    im = im2double(reshape(permute(im, [3 1 2]), 3, []));
end

% Compute the 3D points in world coordinates
[Y, X] = ndgrid(1:size(Z, 1), 1:size(Z, 2));
X = X .* Z;
Y = Y .* Z;
Z = [X(:)'; Y(:)'; Z(:)'];
Z(4,:) = 1;
Z = [P; 0 0 0 1] \ Z;
Z = Z(1:3,:);

% Render
h = render_points(Z, im, num_colors, marker, marker_size);
