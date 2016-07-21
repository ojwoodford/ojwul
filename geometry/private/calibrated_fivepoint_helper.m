%CALIBRATED_FIVEPOINT_HELPER Helper for the calibrated five point algorithm
function varargout = calibrated_fivepoint_helper(varargin)
sourceList = {'calibrated_fivepoint_helper.c'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return