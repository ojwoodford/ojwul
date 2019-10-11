%MESH_AREA Compute the surface area of a triangulated mesh
%
%   [a, A] = mesh_area(V, E)
%
% Compute the surface area of a 3D triangulated mesh
%
%IN:
%   V - 3xM matrix of 3D vertex coordinates.
%   E - 3xN matrix of the 3 vertex indices for each of N triangles.
%
%OUT:
%   a - scalar surface area of entire mesh.
%   A - 1xN surface area of each triangle.

function [a, A] = mesh_area(V, E)
A = normd(cross(V(:,E(2,:)) - V(:,E(1,:)), V(:,E(3,:)) - V(:,E(1,:))), 1);
a = sum(A) * 0.5;
if nargout > 1
    A = A * 0.5;
end
end