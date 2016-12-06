%ADD_DATATIP Add datatip data to graphical objects
%
%   h = add_datatip(h, value)
%   h = add_datatip(h, 'Name', value, ...)
%
% Example:
%   add_datatip(plot(rand(3, 1)), 'Name', {'Tom', 'Dick', 'Harry'});
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
%
%OUT:
%   h - The input scalar graphics handle, for inline usage (see example).

function h = add_datatip(h, varargin)
assert(isscalar(h), 'Only one handle expected');

% Set the 'UserData' property of this object
if nargin == 2
    set(h, 'UserData', varargin{1});
else
    varargin = reshape(varargin, 2, []);
    assert(all(cellfun(@ischar, varargin(1,:))), 'Name, value pairs expected');
    set(h, 'UserData', varargin);
end

% Set the data tip function
set(datacursormode(ancestor(h, 'Figure')), 'UpdateFcn', @datatip_txtfun);
end

function str = datatip_txtfun(~, h)
X = h.Position';
Y = [h.Target.XData; h.Target.YData];
if numel(X) == 3
    Y = [Y; h.Target.ZData];
end
[~, id] = sqdist2closest(X, Y);
data = h.Target.UserData;
str = struct('X', h.Target.XData(id), 'Y', h.Target.YData(id));
if numel(X) == 3
    str.Z = h.Target.ZData(id);
end
if ~isempty(data)
    if iscell(data)
        for a = 1:2:numel(data)
            if iscell(data{a+1})
                str.(data{a}) = data{a+1}{id};
            else
                str.(data{a}) = data{a+1}(:,id)';
            end
        end
    else
        if iscell(data)
            str.Data = data{id};
        else
            str.Data = data(:,id)';
        end
    end
end
str = regexprep(evalc('disp(str)'), '\n *', '\n');
end
