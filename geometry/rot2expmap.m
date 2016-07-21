%ROT2EXPMAP Converts 3x3xN rotation matrices to 3xN exponential mappings
% 
% X = rot2expmap(R)

function X = rot2expmap(R)
sz = [size(R) 1];
th = reshape(acos(1 - mod((0.5 * (1 - R(1,1,:) - R(2,2,:) - R(3,3,:)) - 1), 2)), [1 sz(3:end)]);
X = reshape([R(3,2,:)-R(2,3,:); R(1,3,:)-R(3,1,:); R(2,1,:)-R(1,2,:)], [3 sz(3:end)]);
X = bsxfun(@times, X, th ./ (sqrt(sum(X .* X, 1)) + 1e-300));
