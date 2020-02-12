%SELECT3DPOINTS Select a subset of 3D points
%
%    ind = select3Dpoints(X)
%
% Interactively select a subset of 3D points, by rotating, translating and
% zooming in on the points, then holding down shift and either clicking on
% or drawing round points. Radio buttons allow the user to toggle whether
% selected points are added to or removed from the subset. Close the figure
% to return the selected point indices.
%
%IN:
%   X - 3xN set of points to select from.
%
%OUT:
%   ind - Indices of the selected set of points.

function ind = select3Dpoints(X)
% Initialize variables
global select3DpointsState
select3DpointsState.X = X';
select3DpointsState.M = false(size(X, 2), 1);
select3DpointsState.lastM = select3DpointsState.M;
select3DpointsState.mode = 0;
select3DpointsState.points = zeros(3, 0);

% Create the figure
select3DpointsState.fig = figure('Color', 'w', 'Name', 'Select 3D points by holding down shift and clicking or drawing round points');

% Create the radio buttons
bg = uibuttongroup('Visible', 'off',...
                  'Units', 'pixels', ...
                  'OuterPosition', [0 0 150 70],...
                  'SelectionChangedFcn', @bselection);
uicontrol(bg,'Style',...
             'radiobutton',...
             'String','current OR selection',...
             'UserData', 0, ...
             'Position', [0 50 150 16], ...
             'HandleVisibility','off');
uicontrol(bg,'Style',...
             'radiobutton',...
             'String','current AND selection',...
             'UserData', 1, ...
             'Position', [0 34 150 16], ...
             'HandleVisibility','off');
uicontrol(bg,'Style',...
             'radiobutton',...
             'String','current OR NOT selection',...
             'UserData', 2, ...
             'Position', [0 18 150 16], ...
             'HandleVisibility','off');
uicontrol(bg,'Style',...
             'radiobutton',...
             'String','current AND NOT selection',...
             'UserData', 3, ...
             'Position', [0 2 150 16], ...
             'HandleVisibility','off');
set(bg, 'Visible', 'on');

% Create the plot
select3DpointsState.ax2 = axes('Position', [0 0 1 1]);
select3DpointsState.green = plot(NaN, NaN, 'g-');
xlim([0 1]);
ylim([0 1]);
axis off
select3DpointsState.ax3 = axes('Position', [0 0 1 1]);
select3DpointsState.blue = plot3(select3DpointsState.X(:,1), select3DpointsState.X(:,2), select3DpointsState.X(:,3), 'b.', 'MarkerSize', 1);
hold on;
select3DpointsState.red = plot3(NaN, NaN, NaN, 'r.', 'MarkerSize', 1);
axis equal off
camtarget(median(select3DpointsState.X, 1)');

% Set the callback, pass pointCloud to the callback function
set(select3DpointsState.fig, 'WindowButtonUpFcn', @up, 'WindowButtonDownFcn', @down);
% Require shift to be pressed to select point. Otherwise, control the
% figure.
fcw([0 2 2], '-block');
% Output the selected point indices
ind = find(select3DpointsState.M);
end

function bselection(source, event)
global select3DpointsState
select3DpointsState.mode = event.NewValue.UserData;
end

function render()
global select3DpointsState
X = select3DpointsState.X(~select3DpointsState.M,:);
set(select3DpointsState.blue, 'XData', X(:,1), 'YData', X(:,2), 'ZData', X(:,3));
X = select3DpointsState.X(select3DpointsState.M,:);
set(select3DpointsState.red, 'XData', X(:,1), 'YData', X(:,2), 'ZData', X(:,3));
drawnow();
end

function up(varargin)
global select3DpointsState
set(select3DpointsState.fig, 'WindowButtonMotionFcn', []);
X = select3DpointsState.points';

% Convert the points to 2d
[x, y] = ds2fig(select3DpointsState.ax3, [select3DpointsState.X(:,1); X(:,1)], [select3DpointsState.X(:,2); X(:,2)], [select3DpointsState.X(:,3); X(:,3)]);

% Scale according to dimensions
set(select3DpointsState.ax3, 'Units', 'pixels');
pos = get(select3DpointsState.ax3, 'Position');
x = pos(3) * x;
y = pos(4) * y;
    
n = size(X, 1);
select3DpointsState.lastM = select3DpointsState.M;
if n < 3
    % Select the closest point
    ind = dsearchn([x(1:end-n) y(1:end-n)], [x(end) y(end)]);
    M = false(size(select3DpointsState.M));
    M(ind) = true;
else
    % Select points in the polyshape
    M = inpolygon(x(1:end-n), y(1:end-n), x(end-n+1:end), y(end-n+1:end));
end
% Update the point set
if select3DpointsState.mode > 1
    M = ~M;
end
if mod(select3DpointsState.mode, 2) == 0
    select3DpointsState.M = select3DpointsState.M | M(:);
else
    select3DpointsState.M = select3DpointsState.M & M(:);
end
% Update the rendering
set(select3DpointsState.green, 'XData', NaN, 'YData', NaN);
render();
end

function down(varargin)
global select3DpointsState
X = get(select3DpointsState.ax3, 'CurrentPoint'); 
select3DpointsState.points = X(1,:)';
X = get(select3DpointsState.ax2, 'CurrentPoint'); 
select3DpointsState.selection = X(1,1:2);
set(select3DpointsState.fig, 'WindowButtonMotionFcn', @move);
end

function move(varargin)
global select3DpointsState
X = get(select3DpointsState.ax3, 'CurrentPoint'); 
select3DpointsState.points = [select3DpointsState.points X(1,:)'];
X = get(select3DpointsState.ax2, 'CurrentPoint'); 
select3DpointsState.selection = [select3DpointsState.selection; X(1,1:2)];
set(select3DpointsState.green, 'XData', select3DpointsState.selection(:,1), 'YData', select3DpointsState.selection(:,2));
drawnow();
end