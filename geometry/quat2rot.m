function R = quat2rot(Q)
%QUAT2ROT Convert a quaternion to a 3x3 rotation matrix
%
% R = quat2rot(Q)

R = quat_norm(Q);
R = bsxfun(@times, reshape(R, 4, 1, []), reshape(R, 1, 4, []));
R = [R(1,1,:)+R(2,2,:)-R(3,3,:)-R(4,4,:), 2*(R(2,3,:)-R(4,1,:)), 2*(R(2,4,:)+R(3,1,:));
     2*(R(2,3,:)+R(4,1,:)), R(1,1,:)-R(2,2,:)+R(3,3,:)-R(4,4,:), 2*(R(3,4,:)-R(2,1,:));
     2*(R(2,4,:)-R(3,1,:)), 2*(R(3,4,:)+R(2,1,:)), R(1,1,:)-R(2,2,:)-R(3,3,:)+R(4,4,:)];