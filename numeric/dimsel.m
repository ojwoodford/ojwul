%DIMSEL Select one indexed element from each vector along a given dimension
%
%   B = dimsel(A, I)
%
% Given as input a numeric array, A, and an index array, I, which contains
% indices along a particular dimension of A, output the indexed elements of
% A.
%
% For example, in the code:
%   A = rand(3);
%   [B, I] = max(A, [], 2);
%   C = dimsel(A, I);
% the arrays B and C are identical.
%
%IN:
%   A - Multi-dimensional numeric array.
%   I - Multi-dimensional index array, of the same size as A but for one
%       dimension, which must be singleton. The array must consist of
%       positive intergers in the range [1,N], where N is the length of the
%       dimension of A that is singleton for I.
%
%OUT:
%   B - Output array.

function A = dimsel(A, I)
% Check inputs
szA = size(A);
szI = size(I);
if isequal(szA, szI) || isempty(A)
    return;
end
if numel(szI) ~= numel(szA)
    szI = [szI 1];
end
assert(all(szI == szA | szI == 1), 'Dimensions do not match');
dim = find(szI ~= szA);
assert(numel(dim) == 1, 'Only one dimension of I should be different from that of A');

% Construct the index array
stride = cumprod([1 szA]);
if stride(dim+1) == stride(end)
    I = I(:) * stride(dim) + (1-stride(dim):0)';
elseif stride(dim) == 1
    I = I(:) + (0:stride(dim+1):stride(end)-1)';
else
    I = I(:) * stride(dim) + reshape(bsxfun(@plus, (1-stride(dim):0)', 0:stride(dim+1):stride(end)-1), [], 1);
end

% Construct the output
A = reshape(A(I), szI);
