%INV22N Compute the inverse of an array of 2x2 matrices
%
%   [B, d] = inv22n(A)
%
% Vectorized computation of the inverse of multiple 2x2 matrices.
%
%IN:
%   A - 2x2xN array.
%
%OUT:
%   B - 2x2xN array, where B(:,:,a) = inv(A(:,:,a)).
%   d - 1xN array, where d(a) = det(A(:,:,a)).

function [T, det] = inv22n(T)
sz = size(T);
T = reshape(T, 4, []);
T = T([4 2 3 1],:);
T(2:3,:) = -T(2:3,:);
det = T(1,:) .* T(4,:) - T(2,:) .* T(3,:);
T = reshape(bsxfun(@times, T, 1 ./ det), sz);
end
