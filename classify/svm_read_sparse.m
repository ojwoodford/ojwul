% SVM_READ_SPARSE  Mex wrapper interface to the svm library

function varargout = svm_read_sparse(varargin)
sourceList = {'svm_read_sparse.c'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
