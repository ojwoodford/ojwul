% SVM_PREDICT  Mex wrapper interface to the svm library
% 
%     [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model [,'libsvm_options']);
% 
%         -testing_label_vector:
%             An m by 1 vector of prediction labels. If labels of test
%             data are unknown, simply use any random values.
%         -testing_instance_matrix:
%             An m by n matrix of m testing instances with n features.
%             It can be dense or sparse.
%         -model:
%             The output of svmtrain.
%         -libsvm_options:
%             A string of testing options in the same format as that of LIBSVM.

function varargout = svm_predict(varargin)
sourceDir = 'private/libsvm/';
sourceList = {['-I' sourceDir], 'svm_predict.cpp', [sourceDir 'svm.cpp'], ...
              [sourceDir 'svm_model_matlab.cpp'], '-largeArrayDims'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
