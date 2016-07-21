%KR_from_P Generate K, R and t matrices from a 3x4 projection matrix
%
% [K, R, t] = KR_from_P(P)
%
% P = scalar * K * R * [eye(3) -t]

function [K, R, t] = KR_from_P(P)
st = @(M) M(end:-1:1,end:-1:1)';
[R, K] = qr(st(P(:,1:3)));
K = st(K);
I = diag(K) < 0;
K(:,I) = -K(:,I);
if nargout > 1
    R = st(R);
    R(I,:) = -R(I,:);
    if nargout > 2
        t = (K * R) \ -P(:,4);
    end
end
K = K / K(3,3);