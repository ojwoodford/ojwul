%LINE_PLANE_INTERSECT_SYMEQ Compute the intersection point of a line and a plane
%
%   y = line_plane_intersect_symeq(n, d, x, l)
%
% Symbolically computes the point of intersection of a line and a plane in
% N dimensions.
%
%IN:
%   n - Nx1 plane normal.
%   d - scalar plane offset from origin.
%   x - Nx1 point on line.
%   l - Nx1 line direction.
%
%OUT:
%   y - Nx1 output point.

function y = line_plane_intersect_symeq(n, d, x, l)
persistent eq tv
N = numel(n);
if numel(tv) ~= N * 3 + 1
    % Compute the line-plane intersection equation once
    % Initialize symbolic variables
    n_ = sym(sym('n%d', [N 1]), 'real');
    x_ = sym(sym('x%d', [N 1]), 'real');
    l_ = sym(sym('l%d', [N 1]), 'real');
    syms d_ y_ real
    y_ = solve(n_' * (x_ + l_ * y_) + d_ == 0, y_); % Simultaneous equations
    eq = simplify(x_ + l_ * y_);
    tv = [n_; d_; x_; l_];
end
if nargin < 4
    y = eq;
else
    y = subs(eq, tv, [n; d; x; l]);
end