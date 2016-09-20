%BITWISE_HAMMING Compute all hamming distances between two sets of bit vectors
%
%   C = bitwise_hamming(A, B, [thresh])
%
% Given two sets of bit vectors (each column being a bit vector), compute
% the hamming distances between all pairs of vectors between the two sets.
% If a threshold is given, return only those values below the threshold.
%
%IN:
%   A - MxP matrix of bit vectors.
%   B - MxQ matrix of bit vectors.
%   thresh - A scalar bit count threshold.
%
%OUT:
%   C - PxQ matrix of bitwise hamming distances, or Nx3 is of values be;low
%       the threshold distance, thresh (if given), where each row is 

function C = bitwise_hamming(A, B, thresh)

% Typecast to integer
[A, bytelen] = cast_to_largest_integer(A);
[B, bytelenb] = cast_to_largest_integer(B);
assert(bytelen == bytelenb, 'Input columns must be of the same bit-length');

% Get the input sizes
[a, a] = size(A);
[c, b] = size(B);

% Check if joint array less than 64MB
if a * b * bytelen < 67108864
    % Do the bitxor in one go
    A = repmat(A, [1 1 b]);
    B = repmat(reshape(B, c, 1, b), [1 a 1]);
    C = reshape(bitcount(bitxor(A, B)), a, b);
    if nargin > 2
        M = find(C < thresh);
        L(:,3) = double(C(M));
        [L(:,1), L(:,2)]  = ind2sub(size(C), M);
        C = L;
    end
else
    % Do the bitxor in a loop, to avoid too much memory use
    if nargin > 2
        C = zeros(a*b, 3);
        I = 0;
        if c == 1
            if b > a % Do the smallest for loop
                for i = a:-1:1
                    D = bitcount(bitxor(A(:,i), B))';
                    M = find(D < thresh);
                    if isempty(M)
                        continue;
                    end
                    I = I(end)+1:I(end)+numel(M);
                    C(I,2) = M;
                    C(I,3) = D(M);
                    C(I,1) = i;
                end
            else
                for i = b:-1:1
                    D = bitcount(bitxor(B(:,i), A));
                    M = find(D < thresh);
                    if isempty(M)
                        continue;
                    end
                    I = I(end)+1:I(end)+numel(M);
                    C(I,1) = M;
                    C(I,3) = D(M);
                    C(I,2) = i;
                end
            end
        else
            if b > a % Do the smallest for loop
                for i = a:-1:1
                    D = bitcount(bitxor(repmat(A(:,i), [1 b]), B))';
                    M = find(D < thresh);
                    if isempty(M)
                        continue;
                    end
                    I = I(end)+1:I(end)+numel(M);
                    C(I,2) = M;
                    C(I,3) = D(M);
                    C(I,1) = i;
                end
            else
                for i = b:-1:1
                    D = bitcount(bitxor(repmat(B(:,i), [1 a]), A));
                    M = find(D < thresh);
                    if isempty(M)
                        continue;
                    end
                    I = I(end)+1:I(end)+numel(M);
                    C(I,1) = M;
                    C(I,3) = D(M);
                    C(I,2) = i;
                end
            end
        end
        C = C(1:I(end),:);
    else
        if c == 1
            if b > a % Do the smallest for loop
                for i = a:-1:1
                    C(i,:) = bitcount(bitxor(A(:,i), B));
                end
            else
                for i = b:-1:1
                    C(:,i) = bitcount(bitxor(B(:,i), A))';
                end
            end
        else
            if b > a % Do the smallest for loop
                for i = a:-1:1
                    C(i,:) = bitcount(bitxor(repmat(A(:,i), [1 b]), B));
                end
            else
                for i = b:-1:1
                    C(:,i) = bitcount(bitxor(repmat(B(:,i), [1 a]), A))';
                end
            end
        end
    end
end
end

function [A, bytelen] = cast_to_largest_integer(A)
[a, b] = size(A);
switch class(A)
    case {'double', 'int64', 'uint64'}
        bytelen = 8 * a;
    case {'single', 'int32', 'uint32'}
        bytelen = 4 * a;
    case {'int16', 'uint16'}
        bytelen = 2 * a;
    case {'int8', 'uint8'}
        bytelen = a;
    otherwise
        error('Class %s not supported.', class(A));
end
if mod(bytelen, 8) == 0
    A = reshape(typecast(A(:), 'uint64'), bytelen/8, b);
elseif mod(bytelen, 4) == 0
    A = reshape(typecast(A(:), 'uint32'), bytelen/4, b);
elseif mod(bytelen, 2) == 0
    A = reshape(typecast(A(:), 'uint16'), bytelen/2, b);
else
    A = reshape(typecast(A(:), 'uint8'), bytelen, b);
end
end
