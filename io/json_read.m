function varargout = json_read(varargin)
sourceList = {'json_read.cpp', '-largeArrayDims'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return