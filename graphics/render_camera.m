%RENDER_CAMERA Render 3D camera models, given projection matrices
%
%   h = render_camera(P, [scale, [colour])
%
% Render a set of cameras.
%
%IN:
%   P - 3x4xN camera projection matrices.
%   scale - scalar value indicating the size of the camera. Default: 1.
%   colour - 1x3 colour value for the camera. Default: [0.5 0.5 0.5].
%
%OUT:
%   h - handle to the patch object created.

function h = render_camera(P, scale, colour)

% Create a cone
% Vertices
nSegments = 50;
verts = ones(nSegments+1, 3);
x = linspace(0, pi*2, nSegments+1)';
verts(:,1) = 0.5 * cos(x);
verts(:,2) = 0.5 * sin(x);
verts(end,:) = 0; % Camera center / cone point
verts(:,3) = verts(:,3) * 1.6; % Elongate
% Faces
faces = zeros(nSegments, 3);
faces(:,1) = 1:nSegments;
faces(:,2) = faces(:,1) + 1;
faces(:,3) = nSegments+1;
faces(end,3) = 1;

% Create a cube
% Faces
faces(end+1:end+12,:) = [2 1 3; 2 3 4; 6 5 7; 6 7 8; 5 1 2; 5 2 6; 7 3 4; 7 4 8; 3 1 5; 3 5 7; 4 2 6; 4 6 8] + size(verts, 1);
% Vertices
verts(end+1:end+8,:) = (dec2bin(0:7) - '0') - 0.5;
% Add on a vertical pointer
verts(end+1,:) = [0 -1 0];

% Scale
if nargin > 1
    verts = verts * scale;
end

if nargin > 0 && size(P, 1) == 3 && size(P, 2) == 4
    % Replicate face indices for each camera
    c = numel(P) / 12;
    nv = size(verts, 1);
    faces = reshape(bsxfun(@plus, reshape(faces, [], 1, 3), 0:nv:nv*(c-1)), [], 3);
    % Transform the vertices coordinates into each camera frame
    V = [verts ones(nv, 1)];
    verts = zeros(nv, c, 3);
    for a = 1:c
        % Extract extrinsics
        [K, R, t] = KR_from_P(P(:,:,a));
        % Multiply vertices by inverse
        verts(:,a,:) = reshape(V * [R; t'], [], 1, 3);
    end
    verts = reshape(verts, [], 3);
    
end

% Compute face normals
N = verts';
N = reshape(N(:,faces'), 3, 3, []);
V1 = squeeze(N(:,2,:) - N(:,1,:));
V2 = squeeze(N(:,3,:) - N(:,1,:));
N = V1([2 3 1],:) .* V2([3 1 2],:) - V2([2 3 1],:) .* V1([3 1 2],:);
clear V1 V2
N = bsxfun(@times, N, 1 ./ sqrt(sum(N .* N, 1)));

% Make vertex normals
normals = zeros(size(verts));
normals(faces(:,1),:) = N';

% Render
if nargin < 3
    colour = [0.5 0.5 0.5];
elseif size(colour, 1) == size(P, 3) && size(colour, 2) == 3
    colour = repmat(shiftdim(colour, -1), [nSegments+12 1 1]);
    colour = reshape(colour, [], 3);
end
h = patch('Vertices', verts, 'Faces', faces, 'VertexNormals', normals, 'FaceVertexCData', colour, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BackFaceLighting', 'reverselit', 'FaceLighting', 'flat');
% hold on
% verts = reshape(verts, nSegments+10,[], 3);
% plot3([verts(nSegments+1,:,1); verts(end,:,1)], [verts(nSegments+1,:,2); verts(end,:,2)], [verts(nSegments+1,:,3); verts(end,:,3)], 'r-');