%P_FROM_E Compute potential motion hypotheses from an essential matrix
%
%   P = P_from_E(E)
%
%IN:
%   E - 3x3 essential matrix.
%
%OUT:
%   P - 3x4x4 array of 4 potential extrinsic matrices [R, t] (up to scale).

function P = P_from_E(P)
[U, W, V] = svd(P, 0);
W = [0 -1 0; 1 0 0; 0 0 1];
R = U * W * V';
if det(R) < 0
    U = -U;
    R = -R;
end
U3 = U(:,3);
P = U * W' * V';
P = cat(3, [R U3], [R -U3], [P U3], [P -U3]);
end