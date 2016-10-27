%SKEW Generate 3x3 skew matrices from 3-vectors
%
%   B = skew(A)
%
% Convert one or more 3-vectors to 3x3 skew matrices.
%
%IN:
%   A - 3xM matrix of 3-vectors
%
%OUT:
%   B - 3x3xM array of skew matrices

function B = skew(A)
sz = size(A);
assert(sz(1) == 3);
sz = [3 3 sz(2:end)];
if isa(A, 'sym')
    B = sym(zeros(sz));
else
    B = zeros(sz, class(A));
end
B(2,1,:) = A(3,:);
B(3,1,:) = -A(2,:);
B(1,2,:) = -A(3,:);
B(3,2,:) = A(1,:);
B(1,3,:) = A(2,:);
B(2,3,:) = -A(1,:);
end