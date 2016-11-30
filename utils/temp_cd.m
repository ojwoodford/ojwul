%TEMP_CD Switch to a directory for the duration of the calling function
%
%   temp_cd(dirname)
%
%IN:
%   dirname - Full or relative path to the directory to switch to.

function temp_cd(dirname)
[~, cleanupObj] = fileparts(tempname());
cwd = cd(dirname);
co = onCleanup(@() cd(cwd));
assignin('caller', cleanupObj, co);
end