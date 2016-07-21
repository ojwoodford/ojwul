%EPIPOLAR_DISPARITY Compute the disparity of points on/near epipolar lines
%
%   d = epipolar_disparity(RX, T, x)
%
% Computes the disparity of points on epipolar lines that are closest to
% points in an image.
%
%IN:
%   RX - 3xN world points multiplied by rotation part of projection matrix
%        P(:,1:3). 
%   T - 3x1 or 3xN translation part of projection matrix P(:,4).
%   x - 2xN corresponding points in the image(s).
%
%OUT:
%   d - 1xN output disparity.

function d = epipolar_disparity(RX, T, x)
sum2 = @(X) [X(1,:)+X(2,:); X(3,:)];
RXT = sum2(bsxfun(@times, RX, T));
RX2 = sum2(RX .* RX);
T2 = sum2(T .* T);
Tx = sum(bsxfun(@times, T(1:2,:), x), 1);
RXx = sum(RX(1:2,:) .* x, 1);
d = T(3,:) .* RX2(1,:) - RXT(1,:) .* RX(3,:) + Tx .* RX2(2,:) - RXT(2,:) .* RXx;
d = d ./ (T2(1,:) .* RX(3,:) - RXT(1,:) .* T(3,:) + RXx .* T2(2,:) - RXT(2,:) .* Tx);
end
