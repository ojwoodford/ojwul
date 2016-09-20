% LINEAR_TRAIN  Mex wrapper interface to the linear svm library
%
%  model = linear_train(training_label_vector, training_instance_matrix [,'liblinear_options', 'col']);
% 
%         -training_label_vector:
%             An m by 1 vector of training labels. (type must be double)
%         -training_instance_matrix:
%             An m by n matrix of m training instances with n features.
%             It must be a sparse matrix. (type must be double)
%         -liblinear_options:
%             A string of training options in the same format as that of LIBLINEAR.
%         -col:
%             if 'col' is set, each column of training_instance_matrix is a data instance. Otherwise each row is a data instance.

function varargout = linear_train(varargin)
sourceDir = 'private/liblinear/matlab/';
sourceList = {['-I' sourceDir], 'linear_train.c', [sourceDir 'linear_model_matlab.c'], ...
              [sourceDir '../linear.cpp'], [sourceDir '../tron.cpp'], [sourceDir '../blas/*.c'], ...
              '-largeArrayDims'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
