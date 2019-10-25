%DET33N Compute the determinant of an array of 3x3 matrices
%
%   d = det33n(A)
%
% Vectorized computation of the determinant of multiple 3x3 matrices.
%
%IN:
%   A - 3x3xN array.
%
%OUT:
%   d - 1xN array, where d(a) = det(A(:,:,a)).

function T = det33n(T)
T = reshape(T, 9, []);
cond = zeros(size(T));
d = 1 ./ (1 + abs(T([1 5 9],:)));
cond([1 5 9],:) = d;
T = tmult(reshape(T, 3, 3, []), reshape(cond, 3, 3, []));
T = reshape(T, 9, []);
T = T([1 2 3 1 3 2],:) .* T([5 6 4 6 5 4],:) .* T([9 7 8 8 7 9],:);
T = sum(T(1:3,:) - T(4:6,:));
T = T ./ prod(d);
end
