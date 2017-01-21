%ARRAY_SNAKE_INDICES List of indices snaking through an array
%
%   I = array_snake_indices(sz)
%
% Produces a list of all the indices into an array of size sz, with each
% consecutive index referencing a neighbouring element to the previous
% index.
%
%IN:
%   sz - 1xN vector of the size of array to index.
%
%OUT:
%   I - (prod(sz))x1 uint32 vector of array indices.

function I = array_snake_indices(sz)

% Initialize the indices
I = uint32(1):uint32(prod(sz));
if isempty(I)
    return;
end

% Re-order the indices so they snake through the array
sz = sz(sz>1);
for a = 2:numel(sz)
    I = reshape(I, prod(sz(1:a-1)), sz(a), []);
    I(:,2:2:end,:) = I(end:-1:1,2:2:end,:);
end
I = I(:);

if nargout > 0
    return;
end

% Test the output
X = arrayfun(@(s) {1:s}, sz);
[X{:}] = ndgrid(X{:});
X = cellfun(@(s) col(s(I)), X, 'UniformOutput', false);
X = cat(2, X{:});
Y = sum(abs(diff(X)), 2);
assert(all(Y == 1));
end