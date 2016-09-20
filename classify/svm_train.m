% SVM_TRAIN  Mex wrapper interface to the svm library
%
%   model = svm_train(training_label_vector, training_instance_matrix, [,'libsvm_options']);
% 
%         -training_label_vector:
%             An m by 1 vector of training labels.
%         -training_instance_matrix:
%             An m by n matrix of m training instances with n features.
%             It can be dense or sparse.
%         -libsvm_options:
%             A string of training options in the same format as that of LIBSVM.

function varargout = svm_train(varargin)
sourceDir = 'private/libsvm/';
sourceList = {['-I' sourceDir], 'svm_train.cpp', [sourceDir 'svm.cpp'], ...
              [sourceDir 'svm_model_matlab.cpp'], '-largeArrayDims'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
