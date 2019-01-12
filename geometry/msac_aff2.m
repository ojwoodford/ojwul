% MSAC_AFF2 Fit an affine matrix using RANSAC
%
% Usage:   [A, D] = msac_aff2(x1, x2, t)
%
% Arguments:
%          x1  - 2xN set of homogeneous image points.
%          x2  - 2xN set of image points such that x1<->x2.
%          t   - The distance threshold between data point and the model
%                used to decide whether a point is an inlier or not.
%
% Note that it is assumed that the matching of x1 and x2 are putative and it
% is expected that a percentage of matches will be wrong.
%
% Returns:
%          A - The 3x2 matrix such that A * homg(x1) - x2 is small.
%          D - 1xN array of squared errors.

function [A, D] = msac_aff2(x1, x2, t)
if ~all(size(x1)==size(x2))
    error('Data sets x1 and x2 must have the same dimension');
end

[rows, npts] = size(x1);
if rows ~= 2
    error('x1 and x2 must have 2 rows');
end

% x1 and x2 are 'stacked' to create a 6xN array for ransac
A = msac([x1; x2], @(x) x(3:4,:) / homg(x(1:2,:)), @(A, x) dist(A, x(1:2,:), x(3:4,:)), 3, t);

D = [];
if nargout < 2 || isempty(A)
    return;
end

% Compute the distances
D = dist(A, x1, x2);
end

function D = dist(A, x1, x2)
D = A(:,1:2) * x1 - x2 + A(:,3);
D = sum(D .* D, 1);
end
