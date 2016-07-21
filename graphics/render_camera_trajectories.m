%RENDER_CAMERA_TRAJECTORIES Render the path of a set of cameras
%
%   h = render_camera_tractories(P1, P2, ...)
%   h = render_camera_tractories(P1, I1, ...)
%
% Render a set of camera trajectories
%
%IN:
%   P1 - 3x4xM camera projection matrices defining a camera trajectory.
%   I1 - 1xK indices of keyframes from trajectory P1.
%   P2 - 3x4xN camera projection matrices.
%   etc.
%
%OUT:
%   h - handles to the line objects created for trajectories.

function h = render_camera_trajectories(varargin)

colours = 'brgkmcy';
tf = ishold;
hold on
a = 1;
c = 0;
while a <= nargin
    T = camera_centers(varargin{a});
    h(a) = plot3(T(1,:)', T(2,:)', T(3,:)', [colours(mod(c, numel(colours))+1) '-']);
    a = a + 1;
    if a <= nargin && isvector(varargin{a})
        h(a) = plot3(T(1,varargin{a})', T(2,varargin{a})', T(3,varargin{a})', [colours(mod(c, numel(colours))+1) 'o']);
        a = a + 1;
    end
    c = c + 1;
end
if ~tf
    hold off
end