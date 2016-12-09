%F_FROM_P Compute a fundamental matrix from two projection matrices
%
%   F = F_from_P(P)
%
%IN:
%   P - 3x4x2 array of two projection matrices.
%
%OUT:
%   F - 3x3 fundamental matrix.

function F = F_from_P(F)
X1 = F([2 3],:,1);
X2 = F([3 1],:,1);
X3 = F([1 2],:,1);
Y1 = F([2 3],:,2);
Y2 = F([3 1],:,2);
Y3 = F([1 2],:,2);

F = [det([X1; Y1]) det([X2; Y1]) det([X3; Y1])
    det([X1; Y2]) det([X2; Y2]) det([X3; Y2])
    det([X1; Y3]) det([X2; Y3]) det([X3; Y3])];
end