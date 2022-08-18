%EPIPOLAR_DISPARITY_SYMEQ Compute the disparity of a point on/near an
%                         epipolar line
%
%   d = epipolar_disparity_symeq(RX, T, x)
%
% Symbolically computes the disparity of a point on an epipolar line that
% is closest to a point in an image.
%
%IN:
%   RX - 3x1 world point multiplied by rotation part of projection matrix
%        P(:,1:3). 
%   T - 3x1 translation part of project matrix P(:,4).
%   x - 2x1 corresponding point in the image.
%
%OUT:
%   d - output disparity.

function [d, d_] = epipolar_disparity_symeq(RX, T, x, cov)
persistent eq tv
proj = @(X) X(1:2) / X(3);
if isempty(eq)
    % Compute the point to line distance equation once
    tv = [sym('RX%d', [3 1]); sym('T%d', [3 1]); sym('x%d', [2 1]); col(sym('cov%d', [2 2]))];
    assume(tv, 'real');
    syms d_ real
    eq = (proj(tv(1:3) + tv(4:6) * d_) - tv(7:8));
    eq = simplify(solve(diff(eq' * eq, d_) == sym(0), d_), 'Steps', 10); % Minimize error
end
if nargin < 4
    d = eq;
    if nargout > 1 || nargout == 0
        % Output the point to line residual error
        d_ = simplify((proj(tv(1:3) + tv(4:6) * eq) - tv(7:8)), 'Steps', 10);
        d_ = simplify(sqrt(d_' * d_));
        if nargout == 0
            % Write a MATLAB function to compute the error
            matlabFunction(d_, 'File', 'epipolar_error_symeq');
        end
    end
else
    d = subs(eq, tv, [RX; T; x; col(cov)]);
end
 