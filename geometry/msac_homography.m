% msac_homography - fits fundamental matrix using RANSAC
%
% Usage:   [H, inliers] = msac_homography(x1, x2, t)
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
%          H - The 3x3 homography such that x2 = Hx1.
%          D - 1xN array of squared errors.
%

function [H, D, stats] = msac_homography(x1, x2, t)
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
else
    x1 = hnormalise(x1);
    x2 = hnormalise(x2);
end

% Normalize each set of points so that the origin is at centroid and
% mean distance from origin is sqrt(2).  normalize2dpts also ensures the
% scale parameter is 1.
[T, T] = normalise2dpts([x1 x2]);
x1n = T * x1;
x2n = T * x2;
t = t * T(1) * T(5);

% x1 and x2 are 'stacked' to create a 6xN array for ransac
[H, D, stats] = msac([x1n; x2n], @(x) homography2d(x(1:3,:), x(4:6,:)), @(H, x) homogdist2d(H, x(1:3,:), x(4:6,:)), 4, t);
D = D < t;

% Iterate the least squares fit on the data points considered to
% be inliers
for a = 1:5
    D_ = D;
    H_ = homography2d(x1n(:,D), x2n(:,D));
    if isempty(H_)
        break;
    end
    D = homogdist2d(H_, x1n, x2n) < t;
    if isequal(D_, D) 
        break;
    end
    if sum(D_) > sum(D)
        break;
    end
    H = H_;
end

% Denormalize
H = T \ H * T;
if nargout < 2
    return;
end

% Compute the distances
D = homogdist2d(H, x1, x2);
end

%----------------------------------------------------------------------
% HOMOGRAPHY2D - computes 2D homography
%
% Usage:   H = homography2d(x1, x2)
%
% Arguments:
%          x1  - 3xN set of homogeneous points
%          x2  - 3xN set of homogeneous points such that x1<->x2
%         
%           x  - If a single argument is supplied it is assumed that it
%                is in the form x = [x1; x2]
% Returns:
%          H - the 3x3 homography such that x2 = H*x1
%
% This code follows the normalised direct linear transformation 
% algorithm given by Hartley and Zisserman "Multiple View Geometry in
% Computer Vision" p92.

function H = homography2d(x1, x2)
H = [];
if iscolinear(x1) || iscolinear(x2)
    return;
end
    
% Note that it may have not been possible to normalise
% the points if one was at infinity so the following does not
% assume that scale parameter w = 1.

Npts = size(x1, 2);
A = zeros(3, Npts, 9);
O = [0 0 0];
for n = 1:Npts
X = x1(:,n)';
x = x2(1,n); y = x2(2,n); w = x2(3,n);
A(1,n,:) = [  O  -w*X  y*X];
A(2,n,:) = [ w*X   O  -x*X];
A(3,n,:) = [-y*X  x*X   O ];
end
A = reshape(A, 3*Npts, 9);

[U, D, V] = svd(A, 0); % 'Economy' decomposition for speed

% Extract homography
H = reshape(V(:,9), 3, 3)';
end
    
%----------------------------------------------------------------------
% Function to evaluate the symmetric transfer error of a homography with
% respect to a set of matched points as needed by RANSAC.

function D = homogdist2d(H, x1, x2)
% Calculate, in both directions, the distances   
D = [x1 - hnormalise(H \ x2); x2 - hnormalise(H * x1)];
D = 0.5 * sum(D .* D, 1);
end
    
%----------------------------------------------------------------------
% Test whether any 3 of the 4 points in each set is colinear. 
function tf = iscolinear(x)
try
    tf = any(normd(cross(x(:,[1 1 1 2]) - x(:,[2 2 3 3]), x(:,[1 1 1 2]) - x(:,[3 4 4 4])), 1) < eps);
catch
    tf = false;
end
end    
