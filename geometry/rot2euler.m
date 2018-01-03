%ROT2EULER Converts 3x3xN rotation matrices to 3xN Euler angle representation
% 
% X = rot2euler(R)
%
% The Euler representation is a [roll; pitch; yaw] set of angles (in
% radians), which are applied about the X, Y and Z axes of the frame
% respectively, and in that order.

function X = rot2euler(R)
X = shiftdim(max(min([R(3,3,:) R(3,1,:) R(1,1,:)], 1), -1), 1);
X(2,:) = -asin(X(2,:));
cospitch = cos(X(2,:));
X(1,:) = acos(X(1,:) ./ cospitch);
X(3,:) = acos(X(3,:) ./ cospitch);