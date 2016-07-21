%PROJ2SO3 Project a 3x3 matrix onto the SO(3) manifold
%
%   B = proj2so3(A)
%
% Projects a 3x3 matrix onto SO(3) (the space of rotation matrices) in a
% least squares way, s.t. B = argmin_{B in SO(3)} || B - A ||_F.

function R = proj2so3(R)
[U, D] = eig(R' * R);
[D, I] = sort(diag(D), 'descend');
D = 1 ./ sqrt(D);
D(3) = D(3) * sign(det(R));
U = U(:,I);
R = R * U * diag(D) * U';