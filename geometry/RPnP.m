%RPNP Fast, robust Perspective-n-Points solver
%
%   [P, err] = RPnP(X, x)
%
% Given N>3 3D points and their corresponding 2D projections, compute the
% P = [R t] rigid transform which minimizes the reprojection error.
%
% This is an implementation of the method described in:
%    "A robust O(n) solution to the perspective-n-point problem."
%    S. Li, C. Xu, and M. Xie. IEEE TPAMI, 34(7):1444–1450, 2012.
%
%IN:
%   X - 3xN world coordinates.
%   x - 2xN projected coordinates, i.e. proj(P * [X; ones(1, N)]).
%
%OUT:
%   P - 3x4 [R, t] rigid transform matrix, where R is a rotation matrix and
%       t is a translation matrix.
%   err - Scalar sum of squared reprojection errors.

function [P, minr] = RPnP(X, xx)
P = [];
minr = Inf;

[nd, n] = size(xx);
xv = normalize([xx; ones(1, n)]);

% Select the points with the smallest dot product
x = sum(bsxfun(@times, xv, reshape(xv, nd+1, 1, [])), 1);
[~, x] = min(x(:));
[i1, i2] = ind2sub([n n], x);

% calculating the rotation matrix of $O_aX_aY_aZ_a$.
p2 = X(:,i2);
p0 = (X(:,i1) + p2) / 2;
Rx = x2R(p2 - p0);

% transform the reference points from original object space 
% to the new coordinate frame  $O_aX_aY_aZ_a$.
Xw = Rx * bsxfun(@minus, X, p0);
if any(isnan(Xw(:))) || rank(Xw) < 2
    return;
end

% Dividing the n-point set into (n-2) 3-point subsets
% and setting up the P3P equations
v1 = xv(:,i1);
v2 = xv(:,i2);
cg1 = v1.' * v2;
sg1 = abs(1 - cg1 * cg1);
D1 = Xw(:,i1) - Xw(:,i2);
D1 = D1' * D1;
  
idx = true(1, n);
idx([i1 i2]) = false;
vi = xv(:,idx).';
cg2 = vi * v1;
cg3 = vi * v2;
sg2 = abs(1 - cg2 .* cg2);
vi = Xw(:,idx);
D2 = bsxfun(@minus, Xw(:,i1), vi);
D2 = sum(D2 .* D2)';
D3 = bsxfun(@minus, Xw(:,i2), vi);
D3 = sum(D3 .* D3)';

A1 = D2 ./ D1;
A2 = A1 * sg1 - sg2;
A3 = cg2 .* cg3 - cg1;
A4 = cg1 * cg3 - cg2;
A6 = (D3 - D2 - D1) ./ (2 * D1);
A7 = 1 - cg1 * cg1 - cg2 .* cg2 + cg1 * cg2 .* cg3 + A6 .* sg1;

D4= [A6.^2-A1.*cg3.^2, 2*(A3.*A6-A1.*A4.*cg3),...
     A3.^2+2*A6.*A7-A1.*A4.^2-A2.*cg3.^2,...
     2*(A3.*A7-A2.*A4.*cg3), A7.^2-A2.*A4.^2];

F7= [4*D4(:,1).^2,...
     7*D4(:,2).*D4(:,1),...
     6*D4(:,3).*D4(:,1)+3*D4(:,2).^2,...
     5*D4(:,4).*D4(:,1)+5*D4(:,3).*D4(:,2),...
     4*D4(:,5).*D4(:,1)+4*D4(:,4).*D4(:,2)+2*D4(:,3).^2,...
     3*D4(:,5).*D4(:,2)+3*D4(:,4).*D4(:,3),...
     2*D4(:,5).*D4(:,3)+D4(:,4).^2,...
     D4(:,5).*D4(:,4)];
D7 = sum(F7, 1);

% retriving the local minima of the cost function.
if any(~isfinite(D7))
    return
end
t2s = roots(D7);

maxreal = max(abs(real(t2s)));
t2s(abs(imag(t2s))/maxreal > 0.001)= [];
t2s= real(t2s);

D6 = (7:-1:1) .* D7(1:7);
F6 = polyval(D6,t2s);
t2s(F6 <= 0) = [];

if isempty(t2s)
    return
end

% calculating the camera pose from each local minimum.
m = length(t2s);
for i = 1:m
    % calculating the rotation matrix
    Rx = x2R(v2 * (cg1 + t2s(i)) - v1);
    
    % calculating c, s, tx, ty, tz
    D = zeros(2*n, 6);
    for j = 1:n
        ui= xx(1,j); vi= xx(2,j);
        xi= Xw(1,j); yi= Xw(2,j); zi= Xw(3,j);
        D(2*j-1,:)= [-Rx(2)*yi+ui*(Rx(8)*yi+Rx(9)*zi)-Rx(3)*zi, ...
            -Rx(3)*yi+ui*(Rx(9)*yi-Rx(8)*zi)+Rx(2)*zi, ...
            -1, 0, ui, ui*Rx(7)*xi-Rx(1)*xi];
        D(2*j, :)= [-Rx(5)*yi+vi*(Rx(8)*yi+Rx(9)*zi)-Rx(6)*zi, ...
            -Rx(6)*yi+vi*(Rx(9)*yi-Rx(8)*zi)+Rx(5)*zi, ...
            0, -1, vi, vi*Rx(7)*xi-Rx(4)*xi];
    end
    D = D.' * D;
    [V, ~] = eig(D);
    if V(end,1) == 0
        continue;
    end
    V = V(:,1); V = V / V(end);
    c = V(1); s = V(2); t = V(3:5);
    
    % calculating the camera pose by 3d alignment
    xi= Xw(1,:); yi= Xw(2,:); zi= Xw(3,:);
    XXcs = [Rx(1)*xi+(Rx(2)*c+Rx(3)*s)*yi+(-Rx(2)*s+Rx(3)*c)*zi+t(1);
            Rx(4)*xi+(Rx(5)*c+Rx(6)*s)*yi+(-Rx(5)*s+Rx(6)*c)*zi+t(2);
            Rx(7)*xi+(Rx(8)*c+Rx(9)*s)*yi+(-Rx(8)*s+Rx(9)*c)*zi+t(3)];
    
    XXc = bsxfun(@times, xv, normd(XXcs, 1));
    
    [R, t] = calcampose(XXc,X);
    
    % Calculate the reprojection error
    r = proj(bsxfun(@plus, R * X, t)) - xx;
    r = r(:)' * r(:);
    
    % Keep the smallest reprojection error
    if r > minr
        continue;
    end
    minr = r;
    P = [R t];
end
end


function [R, t] = calcampose(Y,X)

n= length(Y);
K= eye(n)-ones(n,n)/n;

ux= mean(X,2);
uy= mean(Y,2);
sigmx2= mean(sum((X*K).^2));

SXY= Y*K*(X')/n;
[U, D, V]= svd(SXY);
S = eye(3);
if det(SXY) < 0
    S(3,3) = -1;
end

R = U*S*(V');
c2= trace(D*S)/sigmx2;
t= uy-c2*R*ux;

Z = R(:,3);
if norm(xcross(R(:,1), R(:,2)) - Z) > 2e-2
    R(:,3) = -Z;
end
end

function c = xcross(a,b)
c = [a(2)*b(3)-a(3)*b(2);
     a(3)*b(1)-a(1)*b(3);
     a(1)*b(2)-a(2)*b(1)];
end

function R = x2R(x)
if abs(x(2)) < abs(x(3))
    z = xcross(x, [0; 1; 0]);
    y = xcross(z, x);
else
    y = xcross([0; 0; 1], x);
    z = xcross(x, y);
end
R = [x'/norm(x); y'/norm(y); z'/norm(z)];
end
