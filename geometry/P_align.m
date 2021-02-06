%P_ALIGN Align one scene to another
%
%   [Ptgt, Xtgt] = P_align(Pref, Ptgt, Xtgt, [indices])
%
%IN:
%   Pref - 3x4xN camera poses of reference trajectory.
%   Ptgt - 3x4xN camera poses of trajectory to transform.
%   Xtgt - 3xM array of world coordinates to transform.
%   indices - 2x1 indices of camera frames to use for the alignment.
%             Default: [1 2].
%
%OUT:
%   Ptgt - 3x4xN transformed target trajectory
%   Xtgt - 3xM tranformed world coordinates

function [Ptgt, Xtgt] = P_align(Pref, Ptgt, Xtgt, indices)
if nargin < 4
    indices = [1 2];
end
a = indices(1);
b = indices(2);
N = size(Ptgt, 3);

% Adjust the scale
s = norm(cc(Pref(:,:,a)) - cc(Pref(:,:,b))) / norm(cc(Ptgt(:,:,a)) - cc(Ptgt(:,:,b)));
Ptgt(:,4,:) = Ptgt(:,4,:) * s;
Xtgt = Xtgt * s;

% Align the origin
T = [Ptgt(:,:,a); 0 0 0 1] \ [Pref(:,:,a); 0 0 0 1];
for i = 1:N
    Ptgt(:,:,i) = Ptgt(:,:,i) * T;
end
Xtgt = T(1:3,1:3)' * Xtgt + cc(T);
end

function C = cc(P)
C = P(1:3,1:3)' * -P(1:3,4);
end

