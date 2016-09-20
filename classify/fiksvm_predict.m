% FIKSVM_PREDICT  Mex wrapper interface to the svm library
%
%  Usage: [exact_values, pwconst_values, pwlinear_values,[times]] = ...
%         fiksvm_predict(testing_label_vector, testing_instance_matrix, model,'libsvm_options')
%  
%  Output:
%    exact_values    : predictions using binary search
%    pwconst_values  : approximation using piecewise constant function
%    pwlinear_values : approximation using piecewise linear function
%    [times]         : running times 
%  
%  
%  libsvm_options:
%    -b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0);
%    -v verbose flag         : 0 or 1 (default 0);
%    -n number of bins       : [2,...] (default 100);

function varargout = fiksvm_predict(varargin)
sourceDir = 'private/libsvm/';
sourceList = {['-I' sourceDir], 'fiksvm_predict.cpp', [sourceDir 'fiksvm.cpp']}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
