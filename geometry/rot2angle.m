%ROT2ANGLE Given rotation matrices, compute their angle of rotation (in degrees)
%
%   angle = rot2angle(R)
%
%IN:
%   R - 3x{3,4}xN array of rotation or SO(3) matrices.
%
%OUT:
%   angle - 1xN vector of angles, in degrees.

function R = rot2angle(R)
R = reshape(R, 3*size(R, 2), []);
R = sum(R([1 5 9],:), 1);
R = (180 / pi) * acos(max(min(0.5 * (R - 1), 1), -1));
end