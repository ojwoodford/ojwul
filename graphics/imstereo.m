function imstereo(A, B)
if ischar(A)
    A = imread(A);
end
if nargin > 1
    if ischar(B)
        B = imread(B);
    end
    state.sbs = [A B];
else
    state.sbs = A;
    w = floor(size(A, 2) / 2);
    B = A(:,w+1:end,:);
    A = A(:,1:w,:);
end
state.state = 1;
state.offset = 0;

% Create the first image
state.hIm = imdisp(state.sbs);
state.hAx = get(state.hIm, 'Parent');
% Move it to the external monitor
hFig = gcf;
mp = get(0, 'MonitorPositions');
if size(mp, 1) > 1
    p = get(hFig, 'Position');
    p(1:2) = mp(2, 1:2);
    set(hFig, 'Position', p);
end
drawnow;
% Maximize it
set(hFig, 'MenuBar', 'none', 'Color', 'k');
maximize(hFig);
drawnow;

% Prepare resized images
state.maxoffset = 200;
dims = get(hFig, 'Position');
dims = dims(3:4);
w = max(size(A, 2), size(B, 2));
h = size(state.sbs, 1);
f = min((dims - [state.maxoffset 0])./[w h]);
h = floor(h * f);
state.resized = cell(2, 1);
state.resize{1} = imresize(A, [h floor(w*f)], 'bicubic');
state.resize{2} = imresize(B, [h floor(size(B, 2)*f)], 'bicubic');
state.pos = [0 0 max(size(state.resize{1}, 2), size(state.resize{2}, 2))+state.maxoffset size(state.resize{1}, 1)];
state.pos(1:2) = ceil((dims - state.pos(3:4))/2) - 1;

% Set the callback for image navigation, and save the image data in the figure
set(hFig, 'KeyPressFcn', @keypress_callback, 'UserData', state, 'Interruptible', 'off', 'BusyAction', 'cancel');
return

% Keypress callback
% The function which does all the display stuff
function keypress_callback(fig, event_data)
% Check what key was pressed and update the image index as necessary
newstate = 0;
off = 0;
flip = 0;
switch event_data.Character
    case 27 % Escape
        close(fig);
        return
    case 28 % Left
        off = -1;
    case 29 % Right
        off = 1;
    case 30 % Up
        off = -10;
    case 31 % Down
        off = 10;
    case 's'
        newstate = 1;
    case 'h'
        newstate = 2;
    case 'v'
        newstate = 3;
    case 'f'
        flip = 1;
    otherwise
        return
end
% Get the state data
state = get(fig, 'UserData');
if newstate == 0
    newstate = state.state;
else
    state.state = newstate;
end
if flip
    state.resize(1:2) = state.resize([2 1]);
end
state.offset = max(min(state.offset + off, state.maxoffset), -state.maxoffset);
% Display the new image
switch newstate
    case 1 % Side by side
        A = state.sbs;
        set(state.hAx);
        set(state.hAx, 'Units', 'normalized', 'Position', [0 0 1 1]);
    case 2 % Horizontal stripes
        A = zeros(state.pos(4), state.pos(3), size(state.resize{1}, 3), class(state.resize{1}));
        A(1:2:end,(1:size(state.resize{1}, 2))+floor(state.maxoffset/2)-floor(state.offset/2),:) = state.resize{1}(1:2:end,:,:);
        A(2:2:end,(1:size(state.resize{2}, 2))+floor(state.maxoffset/2)+ceil(state.offset/2),:) = state.resize{2}(2:2:end,:,:);
        set(state.hAx, 'Units', 'pixels', 'Position', state.pos);
    case 3 % Vertical subpixel interleaving
        error('Vertical subpixel interleaving not yet supported');
end
% Set the image data
set(state.hIm, 'CData', A);
% Reset the axes limits
set(state.hAx, 'XLim', [0.5 size(A, 2)+0.5], 'YLim', [0.5 size(A, 1)+0.5]);
drawnow;
% Save the current state
set(fig, 'UserData', state);
return