%CAMERA_CENTERS Compute camera centers from projection matrices
%
%   C = camera_extrinsics(P)
%
% Given an array of 3x4 camera projection matrices, this function computes
% their camera centers in world coordinates.
%
%IN:
%   P - 3x4xM array of camera projection matrices.
%
%OUT:
%   C - 3xM array of camera center vectors matrices.

function C = camera_centers(P)
for a = size(P, 3):-1:1
    [K, R, C(:,a)] = KR_from_P(P(:,:,a));
end
