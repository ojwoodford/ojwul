%JSON_READ Mex wrapper to C++ class for fast reading of JSON files
%
%   var = json_read(filename)   
%
% This function wraps the JSON for Modern C++ class in a mex wrapper, for
% fast loading of JSON files into MATLAB variables.
%
%IN:
%   filename - String of filename (if in current directory) or full or
%              relative path to the JSON file to be read in.
%
%OUT:
%   var - MATLAB interpretation of the JSON data.

% Many thanks to Niels Lohmann for his JSON for Modern C++ class:
%    https://nlohmann.github.io/json/
% Thanks also to Sam Hare for pointing me to said class.

function varargout = json_read(varargin)
sourceList = {'json_read.cpp', '-largeArrayDims'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return