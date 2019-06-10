%COMPUTE_HOMOGRAPHY  Compute a best fit homography from correspondences
%
%   H = compute_homography(X, Y)
%
% Compute a hmography from correspondences using the normalized DLT
% algorithm as described in "Multiple View Geometry in Computer Vision".
%
%IN:
%   X - 2xN array of points in image one.
%   Y - 2xN array of points in image two.
%
%OUT:
%   H - 3x3 best fit homography from X to Y.

function H = compute_homography(X, Y)
% Normalize
[X, Tx] = whiten(X);
[Y, Ty] = whiten(Y);

% Create the constraints matrix
O = [0 0 0];
N = size(X, 2);
constraints = zeros(3 * N, 9);
for i = 1:N
    X_ = [X(:,i)' 1];
    x = Y(1,i);
    y = Y(2,i);
    constraints(i*3+(-2:0),:) = [  O  -X_   y*X_;
                                   X_   O  -x*X_;
                                 -y*X_ x*X_   O ];
end

% Decompose
[~, ~, V] = svd(constraints, 0);

% Extract homography
H = reshape(V(:,9), 3, 3)';

% Denormalise
H = Ty \ H * Tx;
end