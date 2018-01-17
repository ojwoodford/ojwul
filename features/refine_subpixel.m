%REFINE_SUBPIXEL N-dimensional sub-pixel refinement
%
%   [offset, val] = refine_subpixel(A, M)
%
% Computes the offsets and values of the refined positions of maxima/minima
% in an N-dimensional array, by fitting a quadratic around the points.
%
%IN:
%   A - An N-dimensional array of size sz.
%   M - A binary mask of size sz indicating which points in the array are
%       to be refined, or a 1xL list of the indices of such points.
%
%OUT:
%   offset - NxL array of offset positions of the estimated true
%            maximum/minimum from the integer position.
%   val - 1xL estimated values at the subpixel maxima/minima.

function [offset, val] = refine_subpixel(V, c)

% Compute the indices
if islogical(c)
    c = find(c);
end
c = reshape(c, 1, []);

% Compute some dimension values
sz = cumprod(size(V));
sz = [1 sz(1:end-1)]';

% Compute the Jacobian (central differences along each dimension)
Vp = V(bsxfun(@plus, c, sz));
Vn = V(bsxfun(@minus, c, sz));
J = 0.5 * (Vp - Vn);

% Compute the Hessian (no skew though, to ensure positive semi-definite
Vc = V(c);
H = Vp + Vn - 2 * Vc;

% Compute the offsets and values
offset = -J ./ H;
if nargout > 1
    val = col(Vc, 2) + dot(J, offset);
end
assert(all(abs(offset(:)) < 1));
end