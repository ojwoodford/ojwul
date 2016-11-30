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

function [J, c] = refine_subpixel(V, c)

% Compute the indices
if islogical(c)
    c = find(c);
end
c = reshape(c, 1, []);

% Compute some dimension values
sz = cumprod(size(V));
sz = [1 sz(1:end-1)]';
sz2 = bsxfun(@plus, sz, sz') - diag(sz);
sz3 = bsxfun(@minus, sz, sz');

% Compute the Jacobian
J = 0.5 * (V(bsxfun(@plus, c, sz)) - V(bsxfun(@minus, c, sz)));

% Compute the Hessian
c = reshape(c, 1, 1, []);
H = V(bsxfun(@plus, c, sz2)) + V(bsxfun(@minus, c, sz2)) - V(bsxfun(@plus, c, sz3)) - V(bsxfun(@minus, c, sz3));
H = bsxfun(@times, H, eye(size(H, 1))*0.75+0.25);

% Compute the offsets and values
if nargout > 1
    for a = 1:numel(c)
        X = -H(:,:,a) \ J(:,a);
        c(a) = V(c(a)) + J(:,a)' * X;
        J(:,a) = X;
    end
    c = reshape(c, 1, []);
else
    for a = 1:numel(c)
        J(:,a) = -H(:,:,a) \ J(:,a);
    end
end 