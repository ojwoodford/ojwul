%JSON_WRITE Write a MATLAB variable to a JSON file
%
%   json_write(var, fname)   
%
% This function wraps the JSON for Modern C++ class in a mex wrapper, for
% fast writing of MATLAB variables into JSON files.
%
%IN:
%   var - MATLAB variable to be written to a JSON file.
%   filename - String of filename (if in current directory) or full or
%              relative path to the JSON file to be written to.

% Many thanks to Niels Lohmann for his JSON for Modern C++ class:
%    https://nlohmann.github.io/json/
% Thanks also to Sam Hare for pointing me to said class.

function varargout = json_write(varargin)
sourceList = {'json_write.cpp', '-largeArrayDims'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return