%RECURSE_SUBDIRS Run a function recursively on a directory structure
%
%   varargout = recurse_subdirs(func, base)
%
% This function calls a function, passing in the path to each subdirectory
% in the tree of the current directory (i.e. including subdirectories of
% subdirectories).
%
%IN:
%   func - A handle for the function to be called recursively. The function
%          should have the form: varargout = func(base).
%   base - The path to the directory the function is to be run on.
%          Default: cd().
%
%OUT:
%   varargout - A cell array of cell arrays of outputs from each call to
%               func.

function varargout = recurse_subdirs(func, base)
if nargin < 2
    base = cd();
end
% Initialize the output
[varargout{1:nargout}] = deal({});
% Run on this directory
try
    if nargout == 0
        func(base);
    else
        [varargout{1:nargout}] = func(base);
        varargout = cellfun(@(v) {{v}}, varargout);
    end
    fprintf('Successfully done: %s\n', base);
catch
end
% Go through sub directories
for d = dir(base)'
    if d.isdir && d.name(1) ~= '.'
        % Do the recursion
        if nargout == 0
            recurse_subdirs(func, sprintf('%s/%s', base, d.name));
        else
            [v{1:nargout}] = recurse_subdirs(func, sprintf('%s/%s', base, d.name));
            % Store the output
            if ~all(cellfun(@isempty, v))
                if all(cellfun(@isempty, varargout))
                    varargout = v;
                else
                    for a = 1:nargout
                        varargout{a} = cat(2, varargout{a}, v{a});
                    end
                end
            end
        end
    end
end
end