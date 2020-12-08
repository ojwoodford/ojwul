%P_FROM_H Compute motion hypotheses from a homography
%
%   P = P_from_H(H)
%
%IN:
%   H - 3x3 homography on calibrated image coordinates, i.e.
%       K^-1 * X2 = H * K^-1 * X1
%
%OUT:
%   P - 3x4xM array of M potential extrinsic matrices [R, t] (up to scale).
%   N - 3xM array of corresponding plane normals 

function [P, N, H] = P_from_H(H)
% Use the method of Faugeras to compute the rotations
% Do the SVD of H
if det(H) < 0
    H = -H;
end
[U, S, V] = svd(H);

% Normalize the homography
S = diag(S);
H = H / S(2);

% Check for no motion case
if S(1) - S(3) < S(2) * 1e-10
    N = [0; 0; 1];
    P = [H zeros(3, 1)];
    return;
end
V = V';
S2 = S .* S;

% Compute first set
ctheta = 1 ./ (S(2) * (S(1) + S(3)));
stheta = sqrt((S2(1) - S2(2)) * (S2(2) - S2(3))) * ctheta;
ctheta = ctheta * (S2(2) + S(1) * S(3));
R1 = U * [ctheta 0 -stheta; 0 1 0; stheta 0 ctheta] * V;
R2 = U * [ctheta 0 stheta; 0 1 0; -stheta 0 ctheta] * V;
P = cat(3, R1, R2);

% Compute second set
ctheta = 1 ./ (S(2) * (S(1) - S(3)));
stheta = sqrt((S2(1) - S2(2)) * (S2(2) - S2(3))) * ctheta;
ctheta = ctheta * (S(1) * S(3) - S2(2));
R1 = U * [ctheta 0 stheta; 0 -1 0; stheta 0 -ctheta] * V;
R2 = U * [ctheta 0 -stheta; 0 -1 0; -stheta 0 -ctheta] * V;
P2 = cat(3, R1, R2);

% Compute the normals and translations
rank1 = cat(3, H - P, -H - P2);
N = sum(rank1, 1);
n = 1 ./ normd(N, 2);
N = N .* n;
T = sum(rank1 .* N, 2);
N = squeeze(N);

% Ensure the first camera is looking at the plane
P = cat(3, P, P2);
M = N(3,:) > 0;
N(:,M) = -N(:,M);
T(:,:,M) = -T(:,:,M);

% Ensure the second camera is in front of the plane
M = dot(reshape(tmult(P, -T, [1 0]), 3, []), N) < 1;
P(:,4,:) = cat(3, T);
P = P(:,:,M);
N = N(:,M);
end
