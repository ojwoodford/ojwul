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

function [P, N] = P_from_H(H)
% Use the method of Faugeras for homography decomposition
% Do the SVD of H
[U, S, V] = svd(H);

% Check all 3 singular values are sufficiently different
S = diag(S);
if any(S([1 2]) ./ S([2 3]) < 1.00001)
    N = zeros(3, 1);
    P = eye(3, 4);
    return;
end

% Compute first set
s = det(U) * det(V);
V = V';
S2 = S .* S;

tp = sqrt((S2([1 2]) - S2([2 3])) ./ (S2(1) - S2(3)));
tp = bsxfun(@times, tp([1 1 2]), [1 1 -1 -1; 0 0 0 0; -1 1 -1 1]);
T = normalize(U * tp * (S(1) - S(3))); 

ctheta = 1 ./ (S(2) * (S(1) + S(3)));
stheta = sqrt((S2(1) - S2(2)) * (S2(2) - S2(3))) * ctheta;
ctheta = ctheta * (S2(2) + S(1) * S(3));
R1 = s * U * [ctheta 0 -stheta; 0 1 0; stheta 0 ctheta] * V;
R2 = s * U * [ctheta 0 stheta; 0 1 0; -stheta 0 ctheta] * V;

P = cat(3, [R1 T(:,1)], [R2 T(:,2)], [R2 T(:,3)], [R1 T(:,4)]);

% Compute second set
tp(3,:) = -tp(3,:);
T = normalize(U * tp * (S(1) + S(3)));

ctheta = 1 ./ (S(2) * (S(1) - S(3)));
stheta = sqrt((S2(1) - S2(2)) * (S2(2) - S2(3))) * ctheta;
ctheta = ctheta * (S(1) * S(3) - S2(2));
R1 = s * U * [ctheta 0 stheta; 0 -1 0; stheta 0 -ctheta] * V;
R2 = s * U * [ctheta 0 -stheta; 0 -1 0; -stheta 0 -ctheta] * V;

P = cat(3, P, [R1 T(:,1)], [R2 T(:,2)], [R2 T(:,3)], [R1 T(:,4)]);
if nargout < 2
    return
end
% Compute the normal
H = H / median(S);
for a = 8:-1:1
    N(:,a) = P(:,4,a) \ (P(:,1:3,a) - H);
end
end
