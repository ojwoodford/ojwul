%TEMP_CD Switch to a directory for the duration of the calling function
%
%   cwd = temp_cd(dirname)
%
%IN:
%   dirname - Full or relative path to the directory to switch to.
%
%OUT:
%   cwd - Path string to current directory.

function cwd = temp_cd(dirname)
assert(evalin('caller', 'exist(''temp_cd_cleanupObj'', ''var'');') == 0, 'temp_cd cannot be called twice in the same function');
cwd = cd(dirname);
co = onCleanup(@() cd(cwd));
assignin('caller', 'temp_cd_cleanupObj', co);
end