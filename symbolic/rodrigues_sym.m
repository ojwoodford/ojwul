%RODRIGUES_SYM Transform angle-axis to rotation matrix via Rodrigues' formula
%
%   R = rodrigues_sym(axis, angle)
%
%IN:
%   axis - 3x1 normalized rotation axis vector.
%   angle - rotation angle in radians.
%
%OUT:
%   R - 3x3 rotation matrix

function R = rodrigues_sym(axis, angle)
R = sin(angle) * skew(axis) + (1 - cos(angle)) * (axis * axis') + cos(angle) * eye(3);
end
