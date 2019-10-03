function index = clickA3DPoint(pointCloud, faces)
%CLICKA3DPOINT
%   H = CLICKA3DPOINT(POINTCLOUD) shows a 3D point cloud and lets the user
%   select points by clicking on them. The selected point is highlighted 
%   and its index in the point cloud will is printed on the screen. 
%   POINTCLOUD should be a 3*N matrix, represending N 3D points. 
%   Handle to the figure is returned.
%
%   other functions required:
%       CALLBACKCLICK3DPOINT  mouse click callback function
%       ROWNORM returns norms of each row of a matrix
%       
%   To test this function ... 
%       pointCloud = rand(3,100)*100;
%       h = clickA3DPoint(pointCloud);
% 
%       now rotate or move the point cloud and try it again.
%       (on the figure View menu, turn the Camera Toolbar on, ...)
%
%   To turn off the callback ...
%       set(h, 'WindowButtonDownFcn',''); 
%
%   by Babak Taati
%   http://rcvlab.ece.queensu.ca/~taatib
%   Robotics and Computer Vision Laboratory (RCVLab)
%   Queen's University
%   May 4, 2005 
%   revised Oct 30, 2007
%   revised May 19, 2009

% Copyright (c) 2009, Babak Taati
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if nargin < 1
    error('Requires one input arguments.')
end

if size(pointCloud, 1)~=3
    error('Input point cloud must be a 3*N matrix.');
end

% show the point cloud
h = figure();
plot3(pointCloud(1,:), pointCloud(2,:), pointCloud(3,:), 'b.');
cameratoolbar('Show'); % show the camera toolbar
hold on; % so we can highlight clicked points without clearing the figure
if nargin > 1
    patch('Vertices', pointCloud', 'Faces', faces', 'FaceVertexCData', repmat([0.7 0.7 0.7], size(pointCloud, 2), 1), 'EdgeColor', 'k', 'LineStyle', '-', 'FaceColor', 'flat');
end

global pointCloudIndex
pointCloudIndex = [];

% set the callback, pass pointCloud to the callback function
set(h, 'WindowButtonDownFcn', {@callbackClickA3DPoint, pointCloud}, 'UserData', []);
% Require shift to be pressed to select point. Otherwise, control the
% figure.
fcw([0 2 2]);

while isempty(pointCloudIndex)
    pause(0.1);
end
index = pointCloudIndex;
close(h);
end

function callbackClickA3DPoint(src, eventData, pointCloud)
% CALLBACKCLICK3DPOINT mouse click callback function for CLICKA3DPOINT
%
%   The transformation between the viewing frame and the point cloud frame
%   is calculated using the camera viewing direction and the 'up' vector.
%   Then, the point cloud is transformed into the viewing frame. Finally,
%   the z coordinate in this frame is ignored and the x and y coordinates
%   of all the points are compared with the mouse click location and the 
%   closest point is selected.
%
%   Babak Taati - May 4, 2005
%   revised Oct 31, 2007
%   revised Jun 3, 2008
%   revised May 19, 2009

ax = gca();
point = get(ax, 'CurrentPoint'); % mouse click position

% Get 2d coordinates in the current view
[x, y] = ds2fig(ax, [pointCloud(1,:) point(1,1)], [pointCloud(2,:) point(1,2)], [pointCloud(3,:) point(1,3)]);

% Scale according to dimensions
set(ax, 'Units', 'pixels');
pos = get(ax, 'Position');
x = pos(3) * x;
y = pos(4) * y;

% find the nearest neighbour to the clicked point 
global pointCloudIndex
pointCloudIndex = dsearchn([x(1:end-1)' y(1:end-1)'], [x(end) y(end)]);
end
