%QUAT_NORM Normalize quaternions

function Q = quat_norm(Q)
Q = reshape(Q, 4, []);
N = sum(Q .* Q);
M = N == 0;
Q(4,M) = 1;
M = ~M & N ~= 1;
if any(M)
    Q(:,M) = bsxfun(@times, Q(:,M), 1./sqrt(N(M)));
end