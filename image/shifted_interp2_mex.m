function varargout = shifted_interp2_mex(varargin)
sourceList = {'shifted_interp2_mex.cpp', '-Xopenmp'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
