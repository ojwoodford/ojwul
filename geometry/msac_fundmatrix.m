% MSAC_FUNDMATRIX - fits fundamental matrix using RANSAC
%
% Usage:   [F, inliers] = msac_fundmatrix(x1, x2, t)
%
% Arguments:
%          x1  - 2xN or 3xN set of homogeneous points.  If the data is
%                2xN it is assumed the homogeneous scale factor is 1.
%          x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
%          t   - The distance threshold between data point and the model
%                used to decide whether a point is an inlier or not. 
%                Note that point coordinates are normalised to that their
%                mean distance from the origin is sqrt(2).  The value of
%                t should be set relative to this, say in the range 
%                0.001 - 0.01  
%
% Note that it is assumed that the matching of x1 and x2 are putative and it
% is expected that a percentage of matches will be wrong.
%
% Returns:
%          F       - The 3x3 fundamental matrix such that x2'Fx1 = 0.
%          inliers - An array of indices of the elements of x1, x2 that were
%                    the inliers for the best model.

% Copyright (c) 2004-2005 Peter Kovesi

function [F, D, stats] = msac_fundmatrix(x1, x2, t)
if ~all(size(x1)==size(x2))
    error('Data sets x1 and x2 must have the same dimension');
end

[rows,npts] = size(x1);
if rows~=2 && rows~=3
    error('x1 and x2 must have 2 or 3 rows');
end

if rows == 2    % Pad data with homogeneous scale factor of 1
    x1 = [x1; ones(1,npts)];
    x2 = [x2; ones(1,npts)];        
end

% Normalize each set of points so that the origin is at centroid and
% mean distance from origin is sqrt(2).  normalize2dpts also ensures the
% scale parameter is 1.
[T, T] = normalise2dpts([x1 x2]);
x1n = T * x1;
x2n = T * x2;
t = t * T(1) * T(5);

% x1 and x2 are 'stacked' to create a 6xN array for ransac
[F, inliers, stats] = msac([x1n; x2n], @(x) fundmatrix(x(1:3,:), x(4:6,:)), @(F, x) funddist(F, x(1:3,:), x(4:6,:)), 8, t);
inliers = inliers < t;

% Now do a final least squares fit on the data points considered to
% be inliers.
F = fundmatrix(x1n(:,inliers), x2n(:,inliers));

% Denormalize
F = T' * F * T;
F = F / F(3,3);
if nargout < 2
    return;
end

% Compute the distances
D = funddist(F, x1, x2);
end

%--------------------------------------------------------------------------
% FUNDMATRIX - computes fundamental matrix from 8 or more points
%
% Function computes the fundamental matrix from 8 or more matching points in
% a stereo pair of images.  The normalised 8 point algorithm given by
% Hartley and Zisserman p265 is used.  To achieve accurate results it is
% recommended that 12 or more points are used
function F = fundmatrix(x1, x2)
% Build the constraint matrix
F = [x2(1,:)'.*x1(1,:)'   x2(1,:)'.*x1(2,:)'  x2(1,:)' ...
     x2(2,:)'.*x1(1,:)'   x2(2,:)'.*x1(2,:)'  x2(2,:)' ...
     x1(1,:)'             x1(2,:)'            ones(size(x1, 2),1) ];       

[U, D, F] = svd(F,0); % Under MATLAB use the economy decomposition

% Extract fundamental matrix from the column of V corresponding to
% smallest singular value.
F = reshape(F(:,9),3,3)';

% Enforce constraint that fundamental matrix has rank 2 by performing
% a svd and then reconstructing with the two largest singular values.
[U, D, F] = svd(F,0);
F = U * diag([D(1,1) D(2,2) 0]) * F';
end    

%--------------------------------------------------------------------------
% Function to evaluate the first order approximation of the geometric error
% (Sampson distance) of the fit of a fundamental matrix with respect to a
% set of matched points as needed by RANSAC.  See: Hartley and Zisserman,
% 'Multiple View Geometry in Computer Vision', page 270.
function D = funddist(F, x1, x2)
Fx = F * x1;     
x2tFx1 = sum(x2 .* Fx, 1);
Fx = [Fx; F' * x2];

% Evaluate distances
D = (x2tFx1 .* x2tFx1) ./ sum(Fx .* Fx, 1);
end
