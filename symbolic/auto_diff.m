%AUTO_DIFF Output C code to compute function values and their Jacobian
%
%   C = auto_diff(funcs, vars, curr_val)
%
% This function differentiates a symbolic set of resdiuals with respect to
% some parameters.
%
%IN:
%   funcs - Mx1 vector of functional expressions.
%   vars - Nx1 vector of variables to differentiate over.
%   curr_val - Nx1 vector of current values for each variable.
%
%OUT:
%   C - text array of C code to compute M function outputs and MxN Jacobian
%       coefficients.

function C = auto_diff(funcs, vars, curr_val)

% Compute the function values
F = subs(funcs, vars, curr_val);

% Compute the Jacobian of the functions
J = jacobian(funcs, vars);

% Linearize the Jacobian w.r.t. the variables
J = subs(J, vars, curr_val);

% Turn into some C code
C = write_C_fragment([F(:)' J(:)']); 
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
str = strrep(sprintf('   %s\n', str{:}), 'A0[0][', 'out[');
% Add float to the start of temporary variables
str = regexprep(str, '   t(\d*) = ', '   float t${$1} = ');
% Replace special strings
% Replace brackets
sub1 = @(n) num2str(str2double(n) - 1);
str = regexprep(str, '_l_(\d+)_r_', '[${sub1($1)}]');
% Put in dots
str = regexprep(str, '_d_', '.');
% Find those outputs which haven't been set
I = regexp(str, 'out\[(\d+)\] = ', 'tokens');
I = find(~ismember((0:numel(symFunc)-1)', str2double([I{:}])'));
% Set the zero outputs
if ~isempty(I)
    str = [str sprintf('   out[%d] = 0;\n', I - 1)];
end
end