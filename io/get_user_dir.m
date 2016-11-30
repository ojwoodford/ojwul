%GET_USER_DIR Get a user/computer-specific directory path
%
%   path_str = get_user_dir(name, check_path, [append])
%
% Ask a user to select a specific directory, and store this path. If a
% valid path already exists, use this.
%
%IN:
%   name - Name of the directory to be found.
%   check_path - Handle to function which takes path_str as input and
%                returns a boolean indicating whether the path is valid.
%   append - String to append to path_str before checking and storing.
%            Default: ''.
%
%OUT:
%   path_str - Path to user specific directory.

function path_str = get_user_dir(name, check_path, append)
if nargin < 3
    append = '';
end

% Check if we already have a valid path
tag = [strrep(lower(name), ' ', '_') '_path'];
path_str = user_string(tag);
if check_path(path_str)
    return;
end

% Ask the user to enter the path
help = sprintf('Please select your %s directory.', name);
while 1
    if ismac() % Is a Mac
        % Give separate warning as the uigetdir dialogue box doesn't have a
        % title
        uiwait(warndlg(help))
    end
    path_str = uigetdir('/', help);
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

