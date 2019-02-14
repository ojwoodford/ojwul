function [D, I] = sqdist2closest(A, B)
%SQDIST2CLOSEST Squared Euclidean distance of closest point from B to A
%
%    [D I] = sqdist2closest(A, B)
%
% IN:
%    A - MxJ matrix of columnwise vectors.
%    B - MxK matrix of columnwise vectors.
%
% OUT:
%    D - 1xJ Squared distance from each vector in A to their respective
%        closest vector in B.
%    I - Indices of closest vectors from B.

% Precompute some things to speed up
D = sum(A .* A, 1);
B2 = sum(B .* B, 1);
if numel(A) < numel(B)
    A = A * 2;
else
    B = B * 2;
end

% Find distance to closest model points
if nargout > 1
    I = zeros(size(D));
    for a = 1:size(A, 2)
        [d, I(a)] = min(B2 - A(:,a)' * B);
        D(a) = D(a) + d;
    end
else
    for a = 1:size(A, 2)
        D(a) = D(a) + min(B2 - A(:,a)' * B);
    end
end