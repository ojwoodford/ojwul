%QDELETE Delete a file if it exists, and no warning if it doesn't
%
%   qdelete(dname)
%
% Quietly delete a file, if it exists. Do not throw a warning if
% it doesn't exist.
%
%IN:
%   fname - String containing the name or path of the file to delete.
%
% See also DELETE.

% Copyright Oliver Woodford 2025

function qdelete(fname)
if exist(fname, 'file')
    delete(fname);
end
end
