%EIG22N Compute the eigenvalues and eigenvectors of 2x2 matrices
%
%   [e, V] = eig22n(A)
%
% Vectorized computation of the eigenvalues and eigenvectors of multiple
% 2x2 matrices.
%
%IN:
%   A - 2x2xN array.
%
%OUT:
%   e - 2x1xN array, where e(:,a) = eig(A(:,:,a)).
%   V - 2x2xN array, where [V(:,:,a), ~] = eig(A(:,:,a)).

function [e, V] = eig22n(A)
% First compute the eigenvalues using the quadratic equation
b = A(1,1,:) + A(2,2,:);
e = A(2,1,:) .* A(1,2,:) - A(1,1,:) .* A(2,2,:);
e = sqrt(b .* b + 4 * e);
e = [b + e; b - e] * 0.5;

if nargout < 2
    return;
end

% Now solve for the eigenvectors
V = normalize(reshape(sum(A, 1), 2, 1, []) - reshape(e, 1, 2, []), 1);
V = [-V(2,:,:); V(1,:,:)];
end
