%LINE_PLANE_INTERSECT Compute the intersection point of a line and a plane
%
%   Z = line_plane_intersect(N, X, Y, D)
%
% Computes the point of intersection of a line and a plane in M dimensions.
%
%IN:
%   N - MxN plane normals.
%   X - MxN points on plane.
%   Y - MxN points on line.
%   D - MxN line directions.
%
%OUT:
%   Z - MxN output points.

function Y = line_plane_intersect(N, X, Y, D)
Y = bsxfun(@plus, Y, bsxfun(@times, D, sum(bsxfun(@times, bsxfun(@minus, X, Y), N), 1) ./ sum(bsxfun(@times, D, N), 1)));
end