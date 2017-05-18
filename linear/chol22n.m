%CHOL22N Compute the Cholesky decomposition of an array of 2x2 matrices
%
%   B = chol22n(A)
%
% Vectorized computation of the Cholesky decomposition of multiple 2x2
% matrices.
%
%IN:
%   A - 2x2xN array.
%
%OUT:
%   B - 2x2xN array, where B(:,:,a) = chol(A(:,:,a), 'lower').

% Formula from here: http://metamerist.blogspot.com/2008/03/googlaziness-cholesky-2x2.html

function T = chol22n(T)
a = sqrt(T(1,1,:));
b = T(2,1,:) ./ a;
c = sqrt(T(2,2,:) - b .* b);
T = [a zeros(1, 1, size(T, 3)); b c];
end
