%QMKDIR Make a directory if it doesn't exist, and no warning if it does
%
%   status = qmkdir(dname)
%
% Quietly make a directory, if it doesn't exist. Do not throw a warning if
% it exists already, or an error if it cannot be created.
%
%IN:
%   dname - String containing the name or path of the directory to make.
%OUT:
%   status - 1 if the directory exists or was created, otherwise zero.
%
% See also MKDIR.

% Copyright Oliver Woodford 2025

function status = qmkdir(dname)
if exist(dname, 'dir')
    status = 1;
    return;
end
status = mkdir(dname);
end
