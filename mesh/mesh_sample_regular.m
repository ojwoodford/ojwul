%MESH_SAMPLE_REGULAR Regularly sample points on a mesh
%
%   [P, I] = mesh_sample_regular(V, E, interval)
%
% Sample points regularly over the surface of a triangulated 3D mesh.
%
%IN:
%   V - 3xM matrix of 3D vertex coordinates.
%   E - 3xN matrix of the 3 vertex indices for each of N triangles.
%   interval - > 0: scalar distance between sample points,
%              < 0: a target number of points.
%
%OUT:
%   P - 3xQ matrix of sampled points.
%   I - 1xQ source face indices.


function [P, J] = mesh_sample_regular(V, E, interval)

J = zeros(1, max(1e7, -interval));
P = zeros(3, numel(J));

if interval < 0
    % Estimate an interval based on target number of points
    interval = sqrt(mesh_area(V, E) / -interval);
end

% For each triangle...
I = 0;
for a = 1:size(E, 2)
    % Get triangle vertices
    V_ = V(:,E(:,a));
    
    % Compute the canonical frame
    C = V_(:,[2 3]) - V_(:,[1 1]);
    N = sum(C .* C, 1);
    N(1) = sqrt(N(1));
    N(2) = sqrt(N(2) - ((C(:,1)' * C(:,2)) / N(1)) .^ 2);
        
    % Compute regular sample points
    S = mod(V(:,1)' * C, interval);
    [S, N] = ndgrid((S(1):interval:N(1))/N(1), (S(2):interval:N(2))/N(2));
    S = [S(:)'; N(:)'];
    S = S(:,sum(S, 1)<1);
    if isempty(S)
        continue;
    end
    
    % Compute the points
    P_ = C(:,1) * S(1,:) + C(:,2) * S(2,:);
        
    % Add the canonical frame origin
    P_ = bsxfun(@plus, P_, V_(:,1));
    
    % Add the points to the list
    I = I(end)+1:I(end)+size(P_, 2);
    P(:,I) = P_;
    J(I) = a;
end
P = P(:,1:I(end));
J = J(1:I(end));
end