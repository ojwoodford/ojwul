% SVM_PRECOMP_MODEL  Mex wrapper interface to the svm library

function varargout = svm_precomp_model(varargin)
sourceDir = 'private/libsvm/';
sourceList = {['-I' sourceDir], 'svm_precomp_model.cpp', [sourceDir 'svm.cpp'], ...
              [sourceDir 'svm_model_matlab.cpp'], [sourceDir 'fiksvm.cpp']}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
