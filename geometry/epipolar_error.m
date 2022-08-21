%EPIPOLAR_ERROR Compute the perpendicular error of points to epipolar lines
%
%   d = epipolar_error(RX, T, x, [c])
%
% Computes the point to line error (i.e. signed distance) of points in an
% image to the corresponding epipolar lines.
%
%IN:
%   RX - 3xN world points multiplied by rotation part of projection matrix
%        P(:,1:3).
%   T - 3x1 or 3xN translation part of projection matrix P(:,4).
%   x - 2xN corresponding points in the image(s).
%   c - 2x2xN conditioning matrix to weight the error.
%
%OUT:
%   d - 1xN output error.

function d = epipolar_error(RX, T, x, cov)
T = proj(T);
x = x - T;
RX = proj(RX) - T;
RX = [RX(2,:); -RX(1,:)];
if nargin > 3
    RX = tmult(inv22n(cov), RX, [1 0]);
    x = tmult(cov, x);
end
d = dot(x, normalize(RX));
end
