function [Da, Ia, Db, Ib] = sqdist2closest(A, B)
%SQDIST2CLOSEST Squared Euclidean distance of closest point from B to A
%
%    [Da, Ia, Db, Ib] = sqdist2closest(A, B)
%
% IN:
%    A - MxJ matrix of columnwise vectors.
%    B - MxK matrix of columnwise vectors.
%
% OUT:
%    Da - 1xJ Squared distance from each vector in A to their respective
%        closest vector in B.
%    Ia - Indices of closest vectors from B.
%    Db - 1xK Squared distance from each vector in B to their respective
%        closest vector in A.
%    Ib - Indices of closest vectors from A.

% Precompute some things to speed up
if numel(A) > numel(B)
    flip = true;
    [A, B] = deal(B, A);
else
    flip = false;
end
Da = sum(A .* A, 1);
A = A * 2;
B2 = sum(B .* B, 1);

Ia = zeros(size(Da));
Db = zeros(1, size(B, 2));
Ib = zeros(size(Db));

% Find distance to closest model points
if flip || nargout > 2
    if nargout > 2
        Db = (B2 - A(:,1)' * B) + Da(1);
        [Da(1), Ia(1)] = min(Db);
        Ib(:) = 1;
        for a = 2:size(A, 2)
            D_ = (B2 - A(:,a)' * B) + Da(a);
            [Da(a), Ia(a)] = min(D_);
            Ib(D_ < Db) = a;
            Db = min(Db, D_);
        end
    elseif nargout > 1
        Db = (B2 - A(:,1)' * B) + Da(1);
        Ib(:) = 1;
        for a = 2:size(A, 2)
            D_ = (B2 - A(:,a)' * B) + Da(a);
            Ib(D_ < Db) = a;
            Db = min(Db, D_);
        end
    else
        Db = (B2 - A(:,1)' * B) + Da(1); 
        for a = 2:size(A, 2)
            Db = min(Db, (B2 - A(:,a)' * B) + Da(a));
        end
    end
elseif nargout > 1
    for a = 1:size(A, 2)
        [d, Ia(a)] = min(B2 - A(:,a)' * B);
        Da(a) = Da(a) + d;
    end
else
    for a = 1:size(A, 2)
        Da(a) = Da(a) + min(B2 - A(:,a)' * B);
    end
end

if flip
    [Da, Ia, Db, Ib] = deal(Db, Ib, Da, Ia);
end
end