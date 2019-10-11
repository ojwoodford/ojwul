%MESH_SAMPLE_RANDOM Randomly sample points on a mesh
%
%   [P, N] = mesh_sample_random(V, E, nPts)
%
% Randomly sample points from a uniform distrubution over the surface area
% of a triangulated 3D mesh.
%
%IN:
%   V - 3xM matrix of 3D vertex coordinates.
%   E - 3xN matrix of the 3 vertex indices for each of N triangles.
%   nPts - integer number of points to sample on the mesh.
%
%OUT:
%   P - 3x(nPts) matrix of sampled points.

function [P, N] = mesh_sample_random(V, E, nPts)

% Compute a cdf of area of triangles
[A, A] = mesh_area(V, E);
A = cumsum(A);

% Sample nPts and assign to triangles
A = histc(rand(nPts, 1), [0 A ./ A(end)]);
A(end-1) = A(end-1) + A(end);
A(end) = 0;

P = zeros(3, nPts);
N = zeros(3, nPts);
nPts = 0;

% For each triangle...
for a = find(A)'    
    % Get triangle vertices
    V_ = V(:,E(:,a));
    
    % Compute the canonical frame
    C = V_(:,[2 3]) - V_(:,[1 1]);
    
    % Compute the normal
    N_ = cross(C(:,1), C(:,2));
    N_ = N_ / norm(N_);
        
    % Compute random sample points
    S = rand(2, A(a));
    M = sum(S, 1) > 1;
    S(:,M) = 1 - S(:,M);
    
    % Compute the points
    P_ = C(:,1) * S(1,:) + C(:,2) * S(2,:);
        
    % Add the canonical frame origin
    P_ = bsxfun(@plus, P_, V_(:,1));
    
    % Add the points to the list
    I = nPts+1:nPts+size(P_, 2);
    P(:,I) = P_;
    N(1,I) = N_(1);
    N(2,I) = N_(2);
    N(3,I) = N_(3);
    nPts = nPts + size(P_, 2);
end
end