%LINE_INTERSECT_SYMEQ Compute the intersection point of two lines
%
%   y = line_intersect_symeq(x1, n1, x2, n2)
%
% Symbolically computes the point of intersection of two lines in 2D.
%
%IN:
%   x1 - 2x1 point on line 1.
%   n1 - 2x1 direction of line 1.
%   x2 - 2x1 point on line 2.
%   n2 - 2x1 direction of line 2.
%
%OUT:
%   v - 2x1 output point.

function y = line_intersect_symeq(x1, n1, x2, n2)
persistent eq tv
if isempty(eq)
    % Compute the line intersection equation once
    syms l1 l2 x1_1 x1_2 n1_1 n1_2 x2_1 x2_2 n2_1 n2_2 real
    tv = [x1_1 x1_2 n1_1 n1_2 x2_1 x2_2 n2_1 n2_2]';
    eq = tv(1:2) + tv(3:4) * l1;
    l1_sol = getfield(solve(eq == tv(5:6) + tv(7:8) * l2, l1, l2), 'l1'); % Simultaneous equations
    %l1_sol = getfield(solve(diff(sum((tv(5:6) + tv(7:8) * l2 - eq) .^ 2), l1) == 0, l1, l2), 'l1'); % Minimize error
    eq = simplify(subs(eq, l1, l1_sol));
end
if nargin < 4
    y = eq;
else
    y = subs(eq, tv, [x1; n1; x2; n2]);
end