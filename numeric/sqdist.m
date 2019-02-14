%SQDIST Squared Euclidean distance between sets of vectors
%
%    D = sqdist(A, [B])
%
% IN:
%    A - MxJ matrix of columnwise vectors.
%    B - MxK matrix of columnwise vectors. Default: B = A.
%
% OUT:
%    D - JxK Squared distance between each vector in A and each vector in B.

function D = sqdist(A, B)
if nargin < 2
    B = A;
    A2 = A .* A;
    B2 = A2;
    A = permute(-2 * A, [2 3 1]);
    A2 = permute(A2, [2 3 1]);
else
    A = permute(A, [2 3 1]);
    A2 = A .* A;
    B2 = B .* B;
    if numel(A) <= numel(B)
        A = -2 * A;
    else
        B = -2 * B;
    end
end
A = reshape(cat(2, ones(size(A)), A, A2), size(A, 1), []);
B = reshape(cat(1, shiftdim(B2, -1), shiftdim(B, -1), ones([1 size(B)])), [], size(B, 2));
D = A * B;
end

