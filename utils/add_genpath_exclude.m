%ADD_GENPATH_EXCLUDE Add a folder and subdirectories to the path, with exclusions
%
%   add_genpath_exclude(folder_path, ...)
%
% For example:
%    add_genpath_exclude('ojwul', '/.git', '\.git')
% adds ojwul and subdirectories to the path, excluding .git folders.
%
%IN:
%   folder_path - Relative or absolute path to the folder to be added to
%                 the path.
%   ... - Each trailing input can be a string to be matched within folders
%         to be excluded.

function add_genpath_exclude(path, varargin)
if ispc()
    token = ';';
else
    token = ':';
end
path = cd(cd(path));
path = genpath(path);
path = strsplit(path, token);
path = path(1:end-1);
for str = varargin
    path = path(cellfun(@isempty, strfind(path, str)));
end
path = sprintf(['%s' token], path{:});
addpath(path);
end
