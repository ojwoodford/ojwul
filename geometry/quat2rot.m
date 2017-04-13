%QUAT2ROT Convert a quaternion to a 3x3 rotation matrix
%
%   R = quat2rot(Q, [w_last])
%
% This function takes in quaternions of the form [w x y z]', such that Q =
% [1 0 0 0]' corresponds to R = eye(3). If w_last is true, quaternions are
% of the form [x y z w].
%
%IN:
%   Q - 4xN array of quaternions in the form [w x y z]'.
%   w_last - boolean indicating whether quaternions should be of the form
%            [x y z w]'. Default: false.
%
%OUT:
%   R - 3x3xN array of corresponding rotation matrices.

function R = quat2rot(Q, w_last)
R = quat_norm(Q);
if nargin > 1 && w_last
    R = R([4 1:3],:);
end
R = bsxfun(@times, reshape(R, 4, 1, []), reshape(R, 1, 4, []));
R = [R(1,1,:)+R(2,2,:)-R(3,3,:)-R(4,4,:), 2*(R(2,3,:)-R(4,1,:)), 2*(R(2,4,:)+R(3,1,:));
     2*(R(2,3,:)+R(4,1,:)), R(1,1,:)-R(2,2,:)+R(3,3,:)-R(4,4,:), 2*(R(3,4,:)-R(2,1,:));
     2*(R(2,4,:)-R(3,1,:)), 2*(R(3,4,:)+R(2,1,:)), R(1,1,:)-R(2,2,:)-R(3,3,:)+R(4,4,:)];
end
