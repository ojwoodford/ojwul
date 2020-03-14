%PROJ2SO3 Project a 3x3 matrix onto the SO(3) manifold
%
%   B = proj2so3(A)
%
% Projects a 3x3 matrix onto SO(3) (the space of rotation matrices) in a
% least squares way, s.t. B = argmin_{B in SO(3)} || B - A ||_F.

function R = proj2so3(R)
R = proj2orthonormal(R);
end