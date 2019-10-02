%GET_USER_PATH Get a user/computer-specific directory or file path
%
%   path_str = get_user_path(name, check_path, type [append])
%
% Ask a user to select a specific directory or file, and store its path.
% If a valid path already exists, use this.
%
%IN:
%   name - Name of the directory to be found.
%   check_path - Handle to function which takes path_str as input and
%                returns a boolean indicating whether the path is valid.
%   type - Scalar. 1 for directory, 2 for file, 3 for application folder.
%   append - String to append to path_str before checking and storing.
%            Default: ''.
%
%OUT:
%   path_str - Path to user specific directory.

function path_str = get_user_path(name, check_path, type, append)
if nargin < 4
    append = '';
end

% Check if we already have a valid path
tag = [strrep(lower(name), ' ', '_') '_path'];
path_str = user_string(tag);
if check_path(path_str)
    return;
end

% Ask the user to enter the path
type_name = {'directory', 'file', 'application folder'};
help = sprintf('Please select your %s %s.', name, type_name{type});
while 1
    if ismac() % Is a Mac
        % Give separate warning as the uigetdir dialogue box doesn't have a
        % title
        uiwait(warndlg(help))
        type = min(type, 2);
    end
    if type == 2
        [file, path_str] = uigetfile('*', help);
        if file == 0
            path_str = 0;
        else
            path_str = [path_str file];
        end
    else
        path_str = uigetdir('/', help);
    end
    if isequal(path_str, 0)
        % User hit cancel or closed window
        error('%s not found.', name);
    end
    path_str = [path_str filesep append];
    
    % Check the path
    if ~check_path(path_str)
        continue;
    end
    
    % Store the valid path
    user_string(tag, path_str);
    return;
end
end

