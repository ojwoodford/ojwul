%COMPUTE_HOMOGRAPHY  Compute a best fit homography from correspondences
%
%   H = compute_homography(X, Y, [normalize])
%
% Compute a hmography from correspondences using the normalized DLT
% algorithm as described in "Multiple View Geometry in Computer Vision".
%
%IN:
%   X - 2xN array of points in image one.
%   Y - 2xN array of points in image two.
%   normalize - boolean flag indicating whether to normalize the inputs.
%
%OUT:
%   H - 3x3 best fit homography from X to Y.

function H = compute_homography(X, Y, normalize)
if nargin < 3
    normalize = true;
end

if normalize
    % Normalize
    [X, Tx] = whiten(X);
    [Y, Ty] = whiten(Y);
end

% Create the constraints matrix
V = zeros(3, size(X, 2), 3, 3);
X = homg(X);
V(:,:,1,2) = -X;
V(:,:,2,1) = X;
Z = X .* Y(1,:);
V(:,:,2,3) = -Z;
V(:,:,3,2) = Z;
Z = X .* Y(2,:);
V(:,:,1,3) = Z;
V(:,:,3,1) = -Z;
V = reshape(permute(V, [3 2 1 4]), [], 9);

% Decompose
[~, ~, V] = svd(V, 0);

% Extract homography
H = reshape(V(:,9), 3, 3)';

if normalize
    % Denormalize
    H = Ty \ H * Tx;
end
end