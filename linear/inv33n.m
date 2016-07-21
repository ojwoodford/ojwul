%INV33N Compute the inverse of an array of 3x3 matrices
%
%   [B, d] = inv33n(A)
%
% Vectorized computation of the inverse of multiple 3x3 matrices.
%
%IN:
%   A - 3x3xN array.
%
%OUT:
%   B - 3x3xN array, where B(:,:,a) = inv(A(:,:,a)).
%   d - 1xN array, where d(a) = det(A(:,:,a)).

function [T, det] = inv33n(T)
sz = size(T);
T = reshape(T, 9, []);
det = T([1 2 3 1 3 2],:) .* T([5 6 4 6 5 4],:) .* T([9 7 8 8 7 9],:);
det = sum(det(1:3,:)) - sum(det(4:6,:));
T = T([5 8 2 7 1 4 4 7 1],:) .* T([9 3 6 6 9 3 8 2 5],:) - T([8 2 5 4 7 1 7 1 4],:) .* T([6 9 3 9 3 6 5 8 2],:);
T = reshape(bsxfun(@times, T, 1 ./ det), sz);
end
