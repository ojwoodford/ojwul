%ADD_DATATIP Add datatip data to graphical objects
%
%   h = add_datatip(h, callback_func, value)
%   h = add_datatip(h, value)
%   h = add_datatip(h, 'Name', value, ...)
%   h = add_datatip(h, callback_func, value)
%   h = add_datatip(h, callback_func, 'Name', value, ...)
%
% Example:
%   add_datatip(plot(rand(3, 1)), @disp, 'Name', {'Tom', 'Dick', 'Harry'});
%
% Add data to graphical objects, for querying using datatips.
%
%IN:
%   h - Scalar graphics object handle, to which to add the data.
%   value - MxN numeric array or 1xN cell array of data, the ith column or
%           cell of which is to be displayed when the ith datapoint is
%           selected.
%   'Name' - String indicating the name of the data value, so that
%            name/value pairs can be used.
%   callback_func - Handle to callback function called with a structure
%                   holding the data as input.
%
%OUT:
%   h - The input scalar graphics handle, for inline usage (see example).

function h = add_datatip(h, varargin)
assert(isscalar(h), 'Only one handle expected');

% Check for a callback
if nargin > 1 && isa(varargin{1}, 'function_handle')
    callback = varargin(1);
    varargin = varargin(2:end);
else
    callback = {};
end

% Check we have name value pairs
if numel(varargin) > 1
    varargin = reshape(varargin, 2, []);
    assert(all(cellfun(@ischar, varargin(1,:))), 'Name, value pairs expected');
end

% Set the 'UserData' property of this object
set(h, 'UserData', [varargin(:); callback], 'Tag', 'add_datatip');

% Set the data tip function
set(datacursormode(ancestor(h, 'Figure')), 'UpdateFcn', @datatip_txtfun);
end

function str = datatip_txtfun(obj, event)
% Get the position of the object that was clicked
X = event.Position';
% Check if the string is cached
persistent last_pos
persistent last_str
if isequal(last_pos, X)
    str = last_str;
    return;
end
% Check if this has an add_datatip structure
id = [];
if strcmp(event.Target.Tag, 'add_datatip')
    % Find the point that was clicked
    hTarg = event.Target;
    if numel(X) == 3
        Y = [hTarg.XData(:)'; hTarg.YData(:)'; hTarg.ZData(:)'];
    else
        Y = [hTarg.XData(:)'; hTarg.YData(:)'];
    end
    [md, id] = sqdist2closest(X, Y);
    Z = Y(:,id);
else
    % Query all the datatip objects
    objs = findobj(ancestor(event.Target, 'Axes'), 'Tag', 'add_datatip');
    Z = NaN(size(X, 1), 1);
    md = Inf;
    for a = 1:numel(objs)
        if numel(X) == 3
            Y = [objs(a).XData(:)'; objs(a).YData(:)'; objs(a).ZData(:)'];
        else
            Y = [objs(a).XData(:)'; objs(a).YData(:)'];
        end
        [d, id] = sqdist2closest(X, Y);
        if d < md
            md = d;
            Z = Y(:,id);
            hTarg = objs(a);
        end
        if md == 0
            break;
        end
    end
end
if isempty(id)
    str = '';
    return;
end
% Check if the string is cached
if isequal(last_pos, Z)
    str = last_str;
    obj.Cursor.Position = Z';
    return;
end
% Generate the structure to display
str = struct('X', Z(1), 'Y', Z(2));
if numel(X) == 3
    str.Z = Z(3);
end
data = hTarg.UserData;
if ~isempty(data)
    try
        N = numel(data) - isa(data{end}, 'function_handle');
        if N > 1
            for a = 1:2:N
                if iscell(data{a+1})
                    str.(data{a}) = data{a+1}{id};
                else
                    str.(data{a}) = data{a+1}(:,id)';
                end
            end
        elseif N == 1
            if iscell(data{1})
                str.Data = data{1}{id};
            else
                str.Data = data{1}(:,id)';
            end
        end
        if N ~= numel(data)
            % Call the callback
            data{end}(str, id);
        end
    catch me
        fprintf('Error in datatip callback: %s\n', getReport(me, 'basic'));
        fprintf('ID: %d\n', id);
        fprintf('Str:\n%s\n', str);
    end
end
% Stringify the structure
str = regexprep(evalc('disp(str)'), '\n *', '\n');
% Store to avoid repetition of computation when we move the cursor
last_pos = Z;
last_str = str;
% Set the data cursor position
obj.Cursor.Position = Z';
end
