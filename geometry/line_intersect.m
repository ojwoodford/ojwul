%LINE_INTERSECT Compute the intersection point of two lines
%
%   [y, l] = line_intersect(x1, n1, x2, n2)
%
% Computes the point of intersection of two lines in 2D.
%
%IN:
%   x1 - 2xN points on line 1.
%   n1 - 2xN directions of line 1.
%   x2 - 2xN points on line 2.
%   n2 - 2xN directions of line 2.
%
%OUT:
%   v - 2xN output point.
%   l - 1xN distance of intersection from x1 in units of norm(n1).

function [y, l] = line_intersect(x1, n1, x2, n2)
n2 = [n2(2,:); -n2(1,:)];
l = (sum(n2 .* x2, 1) - sum(n2 .* x1, 1)) ./ sum(n1 .* n2, 1);
y = x1 + bsxfun(@times, n1, l);