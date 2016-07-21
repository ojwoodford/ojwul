%MAP POINTS  Maps 3D points onto a series of input images
%
%   map_points([options])
%
% Plots the feature points obtained from a Structure from Motion algorithm
% onto the input images in the sequence.
%
% IN:
%   options - Pairs of arguments, the first being the option name and the
%             second being the value. See the m file for details of the
%             options available.

function map_points(varargin)

state.points = 'y.'; % Format string to use when plotting points.
state.tracks = 'c-'; % Format string to use when plotting tracks.
state.depth_thresh = 1e300; % Furthest depth to plot
state.stream_name = 'input.000.png'; % Name of the first frame in the sequence
state.Pi = []; % Camera projection matrices
state.X = []; % 3d scene points
state = vgg_argparse(state, varargin);

% Initialize variables
state.stream = imstream(state.stream_name, 100);
if isempty(state.Pi) || isempty(state.X)
    path = fileparts(state.stream_name);
    if ~isempty(path)
        s = load([path '/data.mat']);
    else
        s = load('data.mat');
    end
    if isempty(state.Pi)
        Pi = PcTo1(s.Pi, state.stream(1)); % Convert P matrices to MATLAB indexing conventions.
    else
        Pi = state.Pi;
    end
    if isempty(state.X)
        points = s.points;
    else
        points = state.X;
    end
else
    Pi = state.Pi;
    points = state.X;
end
state = rmfield(state, {'Pi', 'X'});
if size(points, 1) ~= 3
    points = points';
end
points(4,:) = 1;
state.num_in = size(Pi, 3);
state.ind = 0;
state.projected_points = zeros(state.num_in, size(points, 2), 3);
for n = 1:state.num_in
    % Project the 3-D points
    A = Pi(:,:,n) * points;
    state.projected_points(n,:,3) = A(3,:);
    A(3,:) = 1 ./ A(3,:);
    state.projected_points(n,:,1) = A(1,:) .* A(3,:);
    state.projected_points(n,:,2) = A(2,:) .* A(3,:);
end
if sum(sum(state.projected_points(:,:,3)<0)) > (numel(state.projected_points) / 6)
    state.projected_points(:,:,3) = -state.projected_points(:,:,3);
end
clear A

% Render the figure layout
fig = gcf();
figure(fig); clf('reset');
state.hIm = imdisp(state.stream(1));
hold on
if ~isempty(state.tracks)
    state.hTracks = plot(zeros(100, 1), zeros(100, 1), state.tracks);
end
state.hPoints = plot(zeros(100, 1), zeros(100, 1), state.points);
hold off

% Set the callbacks
set(fig, 'KeyPressFcn', @keypress_callback, 'UserData', state, 'Interruptible', 'off', 'BusyAction', 'cancel');

% Update the figure
event_data.Character = 29;
keypress_callback(fig, event_data);
end

function keypress_callback(fig, event_data)
% Check what key was pressed and update the image index as necessary
switch event_data.Character
    case 28 % Left
        off = -2;
    case 29 % Right
        off = 0;
    case 30 % Up
        off = 9;
    case 31 % Down
        off = -11;
    otherwise
        return
end
state = get(fig, 'UserData');
state.ind = mod(state.ind + off, state.num_in) + 1;

% Select points in front of the camera and within a certain range
A = state.stream(state.ind);
M = state.projected_points(state.ind,:,3) > 0 & state.projected_points(state.ind,:,3) < state.depth_thresh;
M = M & state.projected_points(state.ind,:,1) > 1 & state.projected_points(state.ind,:,1) < size(A, 2);
M = M & state.projected_points(state.ind,:,2) > 1 & state.projected_points(state.ind,:,2) < size(A, 1);

% Update the image and points
set(state.hIm, 'CData', A);
if ~isempty(state.tracks) && state.ind > 1
    X = state.projected_points(max(state.ind-4,1):state.ind,M,1); X(end+1,:) = NaN;
    Y = state.projected_points(max(state.ind-4,1):state.ind,M,2); Y(end+1,:) = NaN;
    set(state.hTracks, 'XData', X(:), 'YData', Y(:));
end
set(state.hPoints, 'XData', state.projected_points(state.ind,M,1), 'YData', state.projected_points(state.ind,M,2));
drawnow;
set(fig, 'UserData', state);
end
   