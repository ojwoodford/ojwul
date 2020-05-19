%FCW  Figure Control Widget: Manipulate figures with key and button presses
%
%   fcw([fig], [buttons], [modifiers])
%   fcw(..., '-link')
%   fcw(..., '-block')
%
% Allows the user to rotate, pan and zoom an axes using key presses and
% mouse gestures. Additionally, press q to quit the widget, r to reset the
% axes and escape to close the figure. This function is non-blocking, but
% fixes axes aspect ratios.
%
% This function plays nicely with other figure callbacks, if called after
% other callbacks have been created. Only figure interactions not used by
% fcw() are passed to previously existing callbacks.
%
% IN:
%   fig - Handle of the figure to be manipulated (default: gcf).
%   buttons - 4x1 cell array indicating the function to associate with
%             each mouse button (left to right) and the scroll action.
%             Functions can be any of:
%                'rot' - Rotate about x and y axes of viewer's coordinate
%                        frame
%                'rotz' - Rotate about z axis of viewer's coordinate frame
%                'zoom' - Zoom (change canera view angle)
%                'zoomz' - Move along z axis of viewer's coordinate frame
%                'pan' - Pan
%                '' - Don't use that button
%             Default: {'rot', 'zoomz', 'pan', 'zoomz'}).
%   modifiers - 1x3 numeric vector indicating which of the modifier keys
%               (shift, ctrl, alt respectively) should be pressed (1), not
%               pressed (0), or can be either (other value), for fcw to
%               react to button/key presses. E.g. [0 1 NaN] indicates that
%               fcw should react to button/key presses if and only if shift
%               is not pressed and ctrl is pressed; alt is ignored.
%               Default: [NaN NaN NaN] (i.e. react to all button presses).
%   '-link' - Specify each axis in the figure to maintain the same pose
%   '-block' - The function doesn't return control to the caller until the
%              widget is quit (by pressing 'q') or the function is closed.
%              This can be useful if collecting data in another callback,
%              to be returned by the calling function.

% (C) Copyright Oliver Woodford 2006-2015

% The initial code came from Torsten Vogel's view3d function, which was in
% turn inspired by rotate3d from The MathWorks, Inc.
% Thanks to Sam Johnson for some bug fixes and good feature requests.

function fcw(varargin)
% Parse input arguments
assert(nargin < 4, 'Too many input arguments');
buttons = {'rot', 'zoomz', 'pan', 'zoom'};
fig = gcf();
link = false;
block = false;
modifiers = NaN(1, 3);
for a = 1:nargin
    v = varargin{a};
    if ischar(v) && numel(v) > 1 && v(1) == '-'
        switch v
            case '-link'
                link = true;
            case '-block'
                block = true;
            otherwise
                error('Unknown option: %s.', v);
        end
    elseif isscalar(v) && ishandle(v)
        fig = ancestor(v, 'figure');
        assert(~isempty(fig), 'Unrecognized handle');
    elseif iscell(v) && numel(v) == 4 && all(cellfun(@ischar, v))
        buttons = v;
    elseif isnumeric(v) && numel(v) == 3
        modifiers = v(:)';
    else
        error('Input not recognized');
    end
end
% Flush any pending draws
drawnow();
% Clear any visualization modes we might be in
pan(fig, 'off');
zoom(fig, 'off');
rotate3d(fig, 'off');
% Find all the 3D axes
hAx = findobj(fig, 'Type', 'axes', '-depth', 1);
M = false(numel(hAx), 1);
for a = 1:numel(hAx)
    zl = zlim(hAx(a));
    M(a) = zl(2) ~= 1 || (zl(1) ~= -1 && zl(1) ~= 0);
end
hAx = hAx(M);
% For each set of axes
data.view = containers.Map('KeyType', 'double', 'ValueType', 'any');
for h = hAx'
    % Set everything to manual
    set(h, 'CameraViewAngleMode', 'manual', 'CameraTargetMode', 'manual', 'CameraPositionMode', 'manual');
    % Store the camera viewpoint
    data.view(double(h)) = camview(h);
end
% Store the data
% Link if necessary
if link
    data.link = linkprop(hAx, {'CameraPosition', 'CameraTarget', 'CameraViewAngle', 'CameraUpVector', 'Projection'});
end
% Get existing callbacks, and quit a previous instance of fcw if running
[prev_mousedown, prev_mouseup, prev_keypress, prev_scroll] = quit_widget(fig);
% Create the data storage object
hData = uipanel(fig, 'Tag', 'fcw_data', 'Visible', 'off', 'Position', [0 0 0 0]);
set(hData, 'UserData', data);
% Create the modifier checking function
M = modifiers == 0 | modifiers == 1;
if ~any(M)
    modifiers = @(varargin) false;
else
    modifiers = modifiers(M);
    mods = {'shift', 'control', 'alt'};
    mods = mods(M);
    modifiers = @(varargin) ~isequal(ismember(mods, get(fig, 'CurrentModifier')), modifiers);
end
% Initialize the callbacks
set(fig, 'WindowButtonDownFcn', [{@fcw_mousedown, {str2func(['fcw_' buttons{1}]), str2func(['fcw_' buttons{2}]), str2func(['fcw_' buttons{3}])}, modifiers} prev_mousedown], ...
         'WindowButtonUpFcn', [{@fcw_mouseup, modifiers} prev_mouseup], ...
         'KeyPressFcn', [{@fcw_keypress, modifiers} prev_keypress], ... 
         'WindowScrollWheelFcn', [{@fcw_scroll, str2func(['fcw_' buttons{4}]), modifiers} prev_scroll], ...
         'BusyAction', 'cancel');
if block
    % Block until the figure is closed or the widget is quit
    while 1
        try
            pause(0.01);
            if isempty(get(hData, 'UserData'))
                break;
            end
        catch
            % Figure was closed
            break;
        end
    end
end
end

function tf = isinvalid(cax, map)
tf = isempty(cax) || ~map.isKey(double(cax)) || ~isequal(size(map(double(cax))), [4 4]);
end

function fcw_keypress(src, eventData, modifiers, varargin)
fig = ancestor(src, 'figure');
cax = get(fig, 'CurrentAxes');
hData = findobj(fig, 'Tag', 'fcw_data', '-depth', 1);
if isinvalid(cax, hData.UserData.view) || modifiers() % Check the required modifiers were pressed, else do nothing
    % Call the other keypress callbacks
    if ~isempty(varargin)
        varargin{1}(src, eventData, varargin{2:end});
    end
    return;
end
step = 1;
if ismember('shift', eventData.Modifier)
    step = 2;
end
if ismember('control', eventData.Modifier)
    step = step * 4;
end
% Which keys do what
switch eventData.Key
    case {'v', 'leftarrow'}
        fcw_pan([], [step 0], cax);
    case {'g', 'rightarrow'}
        fcw_pan([], [-step 0], cax);
    case {'b', 'downarrow'}
        fcw_pan([], [0 step], cax);
    case {'h', 'uparrow'}
        fcw_pan([], [0 -step], cax);
    case {'n', 'x'}
        fcw_rotz([], [0 step], cax);
    case {'j', 's'}
        fcw_rotz([], [0 -step], cax);
    case {'m', 'z'}
        fcw_zoom([], [0 -step], cax);
    case {'k', 'a'}
        fcw_zoom([], [0 step], cax);
    case 'l'
        fcw_rot([], [0 step], cax);
    case 'p'
        fcw_rot([], [0 -step], cax);
    case 'w'
        fcw_rot([], [-step 0], cax);
    case 'q'
        fcw_rot([], [step 0], cax);
    case 'r'
        % Reset all the axes
        for h = findobj(fig, 'Type', 'axes', '-depth', 1)'
            if hData.UserData.view.isKey(double(h))
                camview(h, hData.UserData.view(double(h)));
            end
        end
    case 'q'
        % Quit the widget
        quit_widget(fig);
    case 'escape'
        close(fig);
    otherwise
        % Call the other keypress callbacks
        if ~isempty(varargin)
            varargin{1}(src, eventData, varargin{2:end});
        end
end
end

function fcw_mousedown(src, eventData, funcs, modifiers, varargin)
% Check an axes is selected
fig = ancestor(src, 'figure');
cax = get(fig, 'CurrentAxes');
hData = findobj(fig, 'Tag', 'fcw_data', '-depth', 1);
if isinvalid(cax, hData.UserData.view) || modifiers() % Check the required modifiers were pressed, else do nothing
    % Call the other mousedown callbacks
    % This allows other interactions to be used easily alongside fcw()
    if ~isempty(varargin)
        varargin{1}(src, eventData, varargin{2:end});
    end
    return;
end
% Get the button pressed
switch get(fig, 'SelectionType')
    case 'extend' % Middle button
        method = funcs{2};
    case 'alt' % Right hand button
        method = funcs{3};
    case 'open' % Double click
        camview(cax, hData.UserData.view(double(cax)));
        return;
    otherwise
        method = funcs{1};
end
% Set the cursor
switch func2str(method)
    case {'fcw_zoom', 'fcw_zoomz'}
        shape=[ 2   2   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
                2   1   1   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
                2   1   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
                2   1   2   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
                2   1   2   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN  ;
                2   1   2   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN  ;
                2   1   2   1   1   1   1   1   2 NaN NaN NaN   2   2   2   2  ;
                2   1   2   1   1   2   1   1   1   2 NaN   2   1   2   1   2  ;
                2   1   2   1   2 NaN   2   1   1   1   2   1   1   2   1   2  ;
                2   2   2   2 NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
                NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   2   1   2  ;
                NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
                NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   2   1   2  ;
                NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   1   2  ;
                NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   1   1   2  ;
                NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   2   2  ];
    case 'fcw_pan'
        shape=[ NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
                NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
                NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
                2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
                2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
                NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
                NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
                NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ];
    case {'fcw_rotz', 'fcw_rot'}
        % Rotate
        shape=[ NaN NaN NaN   2   2   2   2   2 NaN   2   2 NaN NaN NaN NaN NaN ;
                NaN NaN NaN   1   1   1   1   1   2   1   1   2 NaN NaN NaN NaN ;
                NaN NaN NaN   2   1   1   1   1   2   1   1   1   2 NaN NaN NaN ;
                NaN NaN   2   1   1   1   1   1   2   2   1   1   1   2 NaN NaN ;
                NaN   2   1   1   1   2   1   1   2 NaN NaN   2   1   1   2 NaN ;
                NaN   2   1   1   2 NaN   2   1   2 NaN NaN   2   1   1   2 NaN ;
                2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
                2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
                2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
                2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
                NaN   2   1   1   2 NaN NaN   2   1   2 NaN   2   1   1   2 NaN ;
                NaN   2   1   1   2 NaN NaN   2   1   1   2   1   1   1   2 NaN ;
                NaN NaN   2   1   1   1   2   2   1   1   1   1   1   2 NaN NaN ;
                NaN NaN NaN   2   1   1   1   2   1   1   1   1   2 NaN NaN NaN ;
                NaN NaN NaN NaN   2   1   1   2   1   1   1   1   1 NaN NaN NaN ;
                NaN NaN NaN NaN NaN   2   2 NaN   2   2   2   2   2 NaN NaN NaN ];
    otherwise
        return
end
% Record where the pointer is
global FCW_POS
FCW_POS = get(0, 'PointerLocation');
% Set the cursor and callback
set(ancestor(src, 'figure'), 'Pointer', 'custom', 'pointershapecdata', shape, 'WindowButtonMotionFcn', {method, cax});
end

function fcw_mouseup(src, eventData, modifiers, varargin)
if modifiers()
    % Call the other mouseup callbacks
    % This allows other interactions to be used easily alongside fcw()
    if ~isempty(varargin)
        varargin{1}(src, eventData, varargin{2:end});
    end
    return;
end
% Clear the cursor and callback
set(ancestor(src, 'figure'), 'WindowButtonMotionFcn', [], 'Pointer', 'arrow');
end

function fcw_scroll(src, eventData, func, modifiers, varargin)
% Get the axes handle
fig = ancestor(src, 'figure');
cax = get(fig, 'CurrentAxes');
hData = findobj(fig, 'Tag', 'fcw_data', '-depth', 1);
if isinvalid(cax, hData.UserData.view) || modifiers() % Check the required modifiers were pressed, else do nothing
    % Call the other mousedown callbacks
    % This allows other interactions to be used easily alongside fcw()
    if ~isempty(varargin)
        varargin{1}(src, eventData, varargin{2:end});
    end
    return;
end
% Call the scroll function
func([], [0 -10*eventData.VerticalScrollCount], cax);
end

function d = check_vals(s, d)
% Check the inputs to the manipulation methods are valid
global FCW_POS
if ~isempty(s)
    % Return the mouse pointers displacement
    new_pt = get(0, 'PointerLocation');
    d = FCW_POS - new_pt;
    FCW_POS = new_pt;
end
end

% Figure manipulation functions
function fcw_rot(s, d, cax)
d = check_vals(s, d);
try
    % Rotate XY
    camorbit(cax, d(1), d(2), 'camera', [0 0 1]);
catch me
    % Error, so release mouse down
    fcw_mouseup(cax);
    fprintf('%s\n', getReport(me, 'extended'));
end
end

function fcw_rotz(s, d, cax)
d = check_vals(s, d);
try
    % Rotate Z
    camroll(cax, d(2));
catch me
    % Error, so release mouse down
    fcw_mouseup(cax);
    fprintf('%s\n', getReport(me, 'extended'));
end
end

function fcw_zoom(s, d, cax)
d = check_vals(s, d);
% Zoom
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2));
try
    camzoom(cax, d);
catch me
    % Error, so release mouse down
    fcw_mouseup(cax);
    fprintf('%s\n', getReport(me, 'extended'));
end
end

function fcw_zoomz(s, d, cax)
d = check_vals(s, d);
% Zoom by moving towards the camera
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2)) - 1;
try
    camdolly(cax, 0, 0, d, 'fixtarget', 'camera');
catch me
    % Error, so release mouse down
    fcw_mouseup(cax);
    fprintf('%s\n', getReport(me, 'extended'));
end
end

function fcw_pan(s, d, cax)
D = check_vals(s, d);
try
    % Pan
    camdolly(cax, D(1), D(2), 0, 'movetarget', 'pixels');
catch me
    % Error, so release mouse down
    fcw_mouseup(cax);
    fprintf('%s\n', getReport(me, 'extended'));
end
end

function A = to_cell(A)
if ~iscell(A)
    if isempty(A)
        A = {};
    else
        A = {A};
    end
end
A = reshape(A, 1, []);
end

function [prev_mousedown, prev_mouseup, prev_keypress, prev_scroll] = quit_widget(fig)
% Get the current callbacks
prev_mousedown = to_cell(get(fig, 'WindowButtonDownFcn'));
prev_mouseup = to_cell(get(fig, 'WindowButtonUpFcn'));
prev_keypress = to_cell(get(fig, 'KeyPressFcn'));
prev_scroll = to_cell(get(fig, 'WindowScrollWheelFcn'));
% Remove the fcw callbacks (assuming they're first)
if ~isempty(prev_mousedown) && strcmp('fcw_mousedown', func2str(prev_mousedown{1}))
    prev_mousedown = prev_mousedown(4:end);
end
if ~isempty(prev_mouseup) && strcmp('fcw_mouseup', func2str(prev_mouseup{1}))
    prev_mouseup = prev_mouseup(2:end);
end
if ~isempty(prev_keypress) && strcmp('fcw_keypress', func2str(prev_keypress{1}))
    prev_keypress = prev_keypress(2:end);
end
if ~isempty(prev_scroll) && strcmp('fcw_scroll', func2str(prev_scroll{1}))
    prev_scroll = prev_scroll(4:end);
end
% Reset the callbacks
set(fig, 'WindowButtonDownFcn', prev_mousedown, 'WindowButtonUpFcn', prev_mouseup, 'KeyPressFcn', prev_keypress, 'WindowScrollWheelFcn', prev_scroll, 'WindowButtonMotionFcn', [], 'Pointer', 'arrow');
% Flag the end
hData = findobj(fig, 'Tag', 'fcw_data', '-depth', 1);
set(hData, 'UserData', []);
delete(hData);
end
