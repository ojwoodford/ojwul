%RENDER_MESH Renders a triangulated mesh
%
%   h = render_mesh(vertices, faces)
%
% Renders a triangulated mesh defined by a set of vertices and faces
% containing vertex indices.
%
%IN:
%   vertices - 
%   faces - 
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

function h = render_mesh(vertices, faces, varargin)
% Set defaults
opts.Marker = 'none';
opts.MarkerSize = 6;
opts.LineStyle = 'none';
opts.FaceColor = 'flat';
opts.LineWidth = 0.5;
% Parse input options
opts = vgg_argparse(opts, varargin);
% Create the options parameters
opts = {'Marker', opts.Marker, 'MarkerSize', opts.MarkerSize, 'LineWidth', opts.LineWidth, 'LineStyle', opts.LineStyle, 'FaceColor', opts.FaceColor};
h = patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', repmat([0.7 0.7 0.7], size(vertices, 1), 1), 'EdgeColor', 'flat', 'MarkerEdgeColor', 'flat', 'BackFaceLighting', 'reverselit', 'FaceLighting', 'phong', opts{:});
