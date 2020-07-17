%ISAPPROX Check if A and B are approximately equal
%
%   [tf, d] = isapprox(A, B, tol)
%
% Determines whether two input arrays are approximately equal. If a
% tolerance is given, and no outputs are requested, the function asserts if
% the inputs aren't approximately equal.
%
%IN:
%   A   - Numeric array.
%   B   - Numeric array to compare against A.
%   tol - Scalar tolerance (proportion of the magnitude).
%
%OUT:
%   tf - Boolean indicating if A and B are approximately equal.
%   d  - Maximum difference between A and B (proportion of the magnitude).
%   D  - Per element different between A and B (proportion of the magnitude).

function [tf, d, D] = isapprox(A, B, tol)
tf = isequal(size(A), size(B));
if ~tf
    if nargout == 0
        error('Matrices have differing sizes: [ %s] and [ %s].', sprintf('%d ', size(A)), sprintf('%d ', size(B)));
    else
        d = -1;
    end
    return;
end
d = abs(A - B);
d2 = abs(A + B);
if issparse(d2)
    [i, j, d] = find(d);
    if isempty(d)
        d = 0;
        d2 = 0;
    else
        i = i + (j - 1) * size(d2, 1);
        d2 = full(d2(i));
    end
end
D = d ./ (d2 + 1);
d = 2 * max(D(:));
if nargin < 3
    if nargout == 0
        fprintf('Matrices differ by scale %g.\n', d);
        return;
    end
    tol = eps * 2;
end
tf = d < tol;
if nargout == 0
    assert(tf, 'Matrices differ outside tolerance. Largest scale difference: %g.', d);
end
end
