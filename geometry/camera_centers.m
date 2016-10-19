%CAMERA_CENTERS Compute camera centers from projection matrices
%
%   C = camera_extrinsics(P, [cam2world])
%
% Given an array of 3x4 camera projection matrices, this function computes
% their camera centers in world coordinates.
%
%IN:
%   P - 3x4xM array of camera projection matrices.
%   cam2world - logical value indicating if projection matrices are camera
%               to world transforms.
%
%OUT:
%   C - 3xM array of camera center vectors matrices.

function C = camera_centers(P, cam2world)
if nargin < 2
    cam2world = false;
end
for a = size(P, 3):-1:1
    [K, R, C(:,a)] = KR_from_P(P(:,:,a));
    if cam2world
        C(:,a) = R * -C(:,a);
    end
end
end
