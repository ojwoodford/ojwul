%EXPMAP2QUAT Convert an exponential mapping rotation to a quaternion
%
% Q = expmap2quat(E)

function Q = expmap2quat(E)

theta = sqrt(sum(E .* E));
M = theta > 1e-8;
m = repmat(0.5, size(theta));
m(M) = sin(theta(M) / 2) ./ theta(M);
Q = [cos(theta / 2); bsxfun(@times, m, reshape(E, 3, []))];