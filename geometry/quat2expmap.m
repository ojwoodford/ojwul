function E = quat2expmap(Q)
%QUAT2EXPMAP Convert a quaternion to an exponential mapping
%
% E = quat2expmap(Q)

theta = acos(Q(1,:));
sintheta = sin(theta);
sintheta(sintheta==0) = 1;
theta = 2 * (theta ./ sintheta);
E = bsxfun(@times, theta, Q(2:4,:));
