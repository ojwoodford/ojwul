%INVSE3N Compute the inverse of an array of 3x4 SE3 matrices
%
%   [B, d] = invSE3n(A)
%
% Vectorized computation of the inverse of multiple 3x4 SE3 matrices.
%
%IN:
%   A - 3x4xN array of SE3 matrices.
%
%OUT:
%   B - 3x4xN array of inverse SE3 matrices.

function T = invSE3n(T)
T(:,1:3,:) = permute(T(:,1:3,:), [2 1 3]); % Transpose the R matrix
T(:,4,:) = tmult(T(:,1:3,:), -T(:,4,:)); % Compute R' * -T
end
