%LINE2LINE_SYMEQ Compute the shortest vector between two lines
%
%   y = line2line_symeq(x1, d1, x2, d2)
%
% Symbolically computes the shortest vector between two lines.
%
%IN:
%   x1 - Nx1 point on line 1.
%   d1 - Nx1 direction vector of line 1.
%   x2 - Nx1 point on line 2.
%   d2 - Nx1 direction vector of line 2.
%
%OUT:
%   y - Nx1 shortest vector.

function y = line2line_symeq(x1, d1, x2, d2)
persistent eq tv
N = numel(x1);
if numel(eq) ~= N
    % Compute the line to line shortest vector once
    % Initialize symbolic variables
    x1_ = sym(sym('x1%d', [N 1]), 'real');
    d1_ = sym(sym('d1%d', [N 1]), 'real');
    x2_ = sym(sym('x2%d', [N 1]), 'real');
    d2_ = sym(sym('d2%d', [N 1]), 'real');
    syms l1 l2 real
    eq = x1_ + d1_ * l1 - x2_ - d2_ * l2;
    l = solve(jacobian(eq' * eq, [l1 l2]) == [0 0], l1, l2); % Simultaneous equations
    eq = simplify(subs(eq, [l1; l2], [l.l1; l.l2]));
    tv = [x1_; d1_; x2_; d2_];
end
if nargin < 4
    y = eq;
else
    y = subs(eq, tv, [x1; d1; x2; d2]);
end