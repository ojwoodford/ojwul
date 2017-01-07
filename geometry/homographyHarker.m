function H = homographyHarker(x1, x2)
%
% Purpose : Computes the general projective transformation between two sets
% of 2D data using a linear algorithm (the homography).
%
% Uses (syntax) :
%   H = homography( DataA, DataB ) 
%
% Input Parameters :
%   DataA, DataB := 2D homogeneous data sets in matrix form (3xn)
%
% Return Parameters :
%   H := the 3x3 homography
%
% Description and algorithms:
%   The algorithm is based on the Direct Linear Transform (DLT) method
%   outlined in Hartley et al.  The method uses orthogonal projections of
%   matrices, such that the vanishing line is treated as the principal
%   component of the reduction.  In this manner, the statistical behaviour
%   of the errors in variables are treated uniformly, see Harker and
%   O'Leary 2005.
%
% References :
%   Harker, M., O'Leary, P., Computation of Homographies, to appear in
%   Proceedings of the British Machine Vision Conference 2005, Oxford,
%   England.
%   Hartley, R., Zisserman, A., Multiple View Geometry in Computer Vision,
%   Cambridge University Press, Cambridge, 2001
%
% Cite this as :
%
% Author : Matthew Harker
% Date : July 25, 2005
% Version : 1.0
%--------------------------------------------------------------------------
% (c) 2005, O'Leary, Harker, University of Leoben, Leoben, Austria
% email: automation@unileoben.ac.at, url: automation.unileoben.ac.at
%--------------------------------------------------------------------------
% History:
%   Date:           Comment:
%   July 25, 2005   Original Version 1.0
%--------------------------------------------------------------------------

% Check input parameters:
[mA, nA] = size(x1);
[mB, nB] = size(x2);
assert(mA == 3 && mB == 3, 'Homogenous coordinates expected');
assert(nA == nB, 'Inputs must have the same size');

% Normalize the input data:
[x1, TA] = normalise2dpts(x1);
[x2, TB] = normalise2dpts(x2);

% Construct the orthogonalized design matrix :
C1 = -x2(1,:) .* x1(1,:);
C2 = -x2(1,:) .* x1(2,:);
C3 = -x2(2,:) .* x1(1,:);
C4 = -x2(2,:) .* x1(2,:);

mC1 = mean(C1);
mC2 = mean(C2);
mC3 = mean(C3);
mC4 = mean(C4);

Mx = [C1' - mC1, C2' - mC2, -x2(1,:)'];
My = [C3' - mC3, C4' - mC4, -x2(2,:)'];

Pp = pinv(x1(1:2,:)');

Bx = Pp * Mx;
By = Pp * My;

D = [Mx - x1(1:2,:)'*Bx; ...
     My - x1(1:2,:)'*By];

% Find v_min and backsubstitute :
[U, S, V] = svd(D, 0);
h789 = V(:,end);
h12 = -Bx * h789;
h45 = -By * h789;
h3 = -[mC1, mC2] * h789(1:2);
h6 = -[mC3, mC4] * h789(1:2);

% Reshape vector h to matrix H, and transform :
H = reshape([h12; h3; h45; h6; h789], 3, 3)';
H = TB \ H * TA;
end
