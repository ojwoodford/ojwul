%EPIPOLAR_ERROR Compute the perpendicular error of points to epipolar lines
%
%   d = epipolar_error(RX, T, x)
%
% Computes the point to line error (i.e. signed distance) of points in an
% image to the corresponding epipolar lines.
%
%IN:
%   RX - 3xN world points multiplied by rotation part of projection matrix
%        P(:,1:3).
%   T - 3x1 or 3xN translation part of projection matrix P(:,4).
%   x - 2xN corresponding points in the image(s).
%
%OUT:
%   d - 1xN output error.

function d = epipolar_error(RX, T, x)
T = proj(T);
RX = normalize(bsxfun(@minus, proj(RX), T));
RX = [RX(2,:); -RX(1,:)];
d = sum(x .* RX, 1) - sum(bsxfun(@times, T, RX), 1);
end
