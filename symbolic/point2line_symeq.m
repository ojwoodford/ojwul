%POINT2LINE_SYMEQ Compute the shortest vector between a point and a line
%
%   V = point2line_symeq(X, Y, D)
%
% Symbolically computes the shortest vector between a point and a line in
% N-D.
%
%IN:
%   X - Nx1 point.
%   Y - Nx1 point on the line.
%   D - Nx1 line direction.
%
%OUT:
%   V - Nx1 output vector.

function v = point2line_symeq(X, Y, D)
persistent eq tv
N = numel(X);
if numel(tv) ~= N * 3
    % Initialize symbolic variables
    tv = sym(sym('tv%d', [N*3 1]), 'real');
    X_ = tv(1:N);
    Y_ = tv(N+1:2*N);
    D_ = tv(2*N+1:3*N);
    % Compute the point to line vector equation once
    syms l real
    eq = X_ - Y_ - D_ * l;
    l_sol = solve(diff(sum(eq .^ 2), l) == sym(0), l); % Minimize error
    eq = simplify(subs(eq, l, l_sol));
end
if nargin < 3
    v = eq;
else
    v = subs(eq, tv, [X; Y; D]);
end