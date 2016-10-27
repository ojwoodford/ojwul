%AUTO_JACOBIAN Write and compile a mex file to compute Jacobian of residuals
%
%   auto_jacobian(residuals, params, fname, var_fixed, var_sum)
%
% This function differentiates a symbolic set of resdiuals with respect to
% some parameters, creates a mex file which can compute the residuals and
% the derivatives, and finally compiles the mex file.
%
%IN:
%   residuals - Mx1 set of symbolic equations for some residuals.
%   params - Nx1 set of symbolic variables to differentiate the residuals
%            with respect to.
%   fname - String containing the file name (or full path) of the mex file
%           to be written, including the extension.
%   var_fixed - Px1 set of symbolic variables which are the same for every
%               column of var_sum variables.
%   var_sum - Qx1 set of symbolic variables, of which many sets can be
%             passed to the mex file (see below).
%
% The resulting mex function operates in the form:
%   [res, J] = fname(var_fixed, var_sum)
%
%IN:
%   var_fixed - Px1 vector of fixed variables.
%   var_sum - QxR array of R vectors of variables which are looped over.
%
%OUT:
%   res - MxR array of residuals
%   J - (M*N)xR array of Jacobians

function auto_jacobian(residuals, params, fname, var_fixed, var_sum)

fprintf('AUTO JACOBIAN: %s\n', fname)

% Residuals with zero inputs
fprintf('Residuals at zero...'); drawnow;
tic
res = subs(residuals, params, zeros(size(params)));
fprintf(' %gs.\n', toc);

% Compute the Jacobian of the residuals
fprintf('Compute Jacobian...'); drawnow;
tic
J = jacobian(residuals, params);
fprintf(' %gs.\n', toc);

% Linearize the Jacobian w.r.t. the parameters
fprintf('Taylor expansion...'); drawnow;
tic
J = subs(J, params, zeros(size(params)));
fprintf(' %gs.\n', toc);

% Turn into some C code
fprintf('Generating C code...'); drawnow;
tic
res = res(:)';
J = J(:)';
res_ = write_C_fragment(res); 
J_ = write_C_fragment(J);

% % Generate the other lines of the C file
lines = { ...
'#include "mex.h"', ...
'#include <math.h>', ...
'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {', ...
'if (nrhs != 2 || nlhs < 1 || nlhs > 2) { mexErrMsgTxt("Unexpected number of arguments."); }', ...
['if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0]) != ' num2str(numel(var_fixed)) ')'], ...
'{ mexErrMsgTxt("First input of unexpected size or type"); }', ...
['if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != ' num2str(numel(var_sum)) ')'], ...
'{ mexErrMsgTxt("Second input of unexpected size or type"); }', ...
'const double *var_fixed = (const double *)mxGetData(prhs[0]);', ...
'const double *var_sum = (const double *)mxGetData(prhs[1]);', ...
'int N = mxGetN(prhs[1]);', ...
['plhs[0] = mxCreateNumericMatrix(' num2str(numel(res)) ', N, mxDOUBLE_CLASS, mxREAL);'], ...
'double *res = (double *)mxGetData(plhs[0]);', ...
'double *J = NULL;', ...
['if (nlhs > 1) { plhs[1] = mxCreateNumericMatrix(' num2str(numel(J)) ', N, mxDOUBLE_CLASS, mxREAL); J = (double *)mxGetData(plhs[1]); }'], ...
['for (int a = 0; a < N; ++a, var_sum += ' num2str(numel(var_sum)) ', res += ' num2str(numel(res)) ') {'], ...
res_, ...
'if (J) { ', ...
J_, ...
['J += ' num2str(numel(J)) ';'], ...
'}}}', ...
};

% Write out to a file
fid = fopen(fname, 'wt');
for l = lines
    fprintf(fid, [l{1} '\n']);
end
fclose(fid);
fprintf(' %gs.\n', toc);

% Mex the function
fprintf('Compiling...'); drawnow;
tic
try
    mex(fname, '-O');
catch me
    % The file might be in use - just report the error
    fprintf('%s\n', getReport(me, 'basic'));
end
fprintf(' %gs.\n', toc);
end

function str = write_C_fragment(symFunc)
% Save the C code to file
fname = [tempname '.cpp'];
ccode(symFunc, 'file', fname);
% Read in the code fragment
fid = fopen(fname, 'rt');
str = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
% Delete the file
delete(fname);
% Deal with single value cases
str = str{1};
if numel(symFunc) == 1
    str{end} = regexprep(str{end}, '\S+ = ', 'A0[0][0] = ');
end
% Remove the extra dimension from the output, and change the variable name
str = strrep(sprintf('   %s\n', str{:}), 'A0[0][', [inputname(1) '[']);
% Add double to the start of temporary variables
str = regexprep(str, '   t(\d*) = ', '   double t${$1} = ');
% Put in the brackets
sub1 = @(n) num2str(str2double(n) - 1);
str = regexprep(str, '_l_(\d*)_r_', '[${sub1($1)}]');
end