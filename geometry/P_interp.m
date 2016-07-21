function P = P_interp(P0, P1, frames)
%P_INTERP  Constant velocity interpolation between 2 projection matrices
%
%   P = P_interp(P0, P1, frames)
%
% Creates a series of projection matrices interpolated from two projection
% matrices. The geodesic path of the cameras is such:
%   - Camera speed (i.e. magnitude of vector velocity) is constant.
%   - Camera angular velocity is constant.
%   - The rotation between the camera frame and local direction of travel
%     is constant.
%
%IN:
%   P0 - 3x4 projection matrix, the first between which the output matrices
%        are interpolated from. 
%   P1 - 3x4 projection matrix, the second between which the output
%        matrices are interpolated from.
%   frames - Nx1 list of frame positions to compute along the geodesic path
%            between P0 and P1, where 0 is the pose of P0 and 1 the pose of
%            P1.
%
%OUT:
%   P - 3x4xN array of output projection matrices.

% 29/05/2014 - Use geodesic path, by doing interpolation on the Lie vector.

% Decompose into constituent parts
[k0, r0, t0] = KR_from_P(P0);
[k1, r1, t1] = KR_from_P(P1);

% If K matrices are identical but for sign, put the sign change into the R
% matrices
a = diag(k0) ./ diag(k1);
b = sign(a);
if all(abs(a - b) < 1e-8) && any(b == -1)
    a = diag(sign(b + 0.5));
    k0 = k0 * a;
    r0 = a * r0;
end
clear a b

% Calculate our step sizes
vec = @(M) [0.5 * (M([7 9 2]) - M([10 3 5]))'; M(1:3,4)];
coord = @(P) real(vec(logm(P)));
P0 = [r0 r0*-t0; 0 0 0 1];
P1 = coord([r1 r1*-t1; 0 0 0 1] / P0);
k1 = k1 - k0;

% Compute the interpolated extrinsics
P = tmult(expm_srt_3d(P1 * frames(:)'), P0);
% Compose with the calibration
P = tmult(bsxfun(@plus, k0, bsxfun(@times, k1, reshape(frames, 1, 1, []))), P);
