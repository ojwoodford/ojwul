% LINEAR_PREDICT  Mex wrapper interface to the linear svm library
%
%  [predicted_label, accuracy, decision_values/prob_estimates] = linear_predict(testing_label_vector, testing_instance_matrix, model [, 'liblinear_options', 'col']);
% 
%         -testing_label_vector:
%             An m by 1 vector of prediction labels. If labels of test
%             data are unknown, simply use any random values. (type must be double)
%         -testing_instance_matrix:
%             An m by n matrix of m testing instances with n features.
%             It must be a sparse matrix. (type must be double)
%         -model:
%             The output of train.
%         -liblinear_options:
%             A string of testing options in the same format as that of LIBLINEAR.
%         -col:
%             if 'col' is set, each column of testing_instance_matrix is a data instance. Otherwise each row is a data instance.

function varargout = linear_predict(varargin)
sourceDir = 'private/liblinear/matlab/';
sourceList = {['-I' sourceDir], 'linear_predict.c', [sourceDir 'linear_model_matlab.c'], ...
              [sourceDir '../linear.cpp'], [sourceDir '../tron.cpp'], [sourceDir '../blas/*.c'], ...
              '-largeArrayDims'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
