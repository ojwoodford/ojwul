%CAMERA_EXTRINSICS Compute camera extrinsics from projection matrices
%
%   [Pe, K] = camera_extrinsics(P)
%
% Given an array of 3x4 camera projection matrices, this function computes
% their extrinsics (matrices in SE3), and (optionally) their intrinsics
% (upper triangular 3x3 matrices).
%
% The output is such that P(:,:,i) = K(:,:,i) * Pe(:,:,i), for all i.
%
%IN:
%   P - 3x4xM array of camera projection matrices.
%
%OUT:
%   Pe - 3x4xM array of SE3 camera extrinsics matrices.
%   K - 3x3xM array of upper triangular camera intrinsics matrices.

function [P, K] = camera_extrinsics(P)
sz = [size(P) 1];
K = zeros([3, 3, sz(3:end)]);
for a = 1:prod(sz(3:end))
    [K(:,:,a), R, t] = KR_from_P(P(:,:,a));
    P(:,:,a) = [R R*-t];
end
