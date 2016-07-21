function Q = quat_mult(Q1, Q2)
%QUAT_MULT Multiply/compose two quaternion rotations or arrays of rotations
%
% Q = quat_mult(Q1, Q2)

Q = bsxfun(@times, reshape(Q1, 4, 1, []), reshape(Q2, 1, 4, []));
Q = squeeze([Q(1,1,:) - Q(2,2,:) - Q(3,3,:) - Q(4,4,:);...
             Q(1,2,:) + Q(2,1,:) + Q(3,4,:) - Q(4,3,:); ...
             Q(1,3,:) - Q(2,4,:) + Q(3,1,:) + Q(4,2,:); ...
             Q(1,4,:) + Q(2,3,:) - Q(3,2,:) + Q(4,1,:)]);
