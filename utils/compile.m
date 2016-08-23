%COMPILE Mex compilation helper function
%
% Examples:
%   compile func1 func2 ... -option1 -option2 ...
%
% This function can be used to (re)compile a number of mex functions, but
% is also a helper function enabling inline compilation.

function varargout = compile(varargin)
% There are two types of call:
% 1. compile('func1', 'func2', ..., [options]).
% 2. [varargout{1:nargout}] = compile(varargin), called from function to
% be compiled.
% Work out which this is

% Try to get the source list
try
    sourceList = evalin('caller', 'sourceList');
catch
    sourceList = [];
end

if iscell(sourceList)
%OJW_MEXCOMPILE_SCRIPT  Compilation helper script for mex files
%
% Should be placed in the m-file of a mexed function, after the
% following lines of code, thus:
%
% function varargout = mymex(varargin)
% sourceList = {'-Iinclude', '-Llib', 'mymex.c', 'mymexhelper.c', '-lmysharedlib'};
% [varargout{1:nargout} = ojw_mexcompile_script(varargin{:});
%
% The script will compile the source inline (i.e. compile & then run) if
% the function has been called without first being compiled.
%
% If varargin{1} == 'compile' and varargin{2} = last_compilation_time,
% then the script calls ojw_mexcompile to compile the function if any of
% the source files have changed since last_compilation_time.
        
    % Get the name of the calling function
    funcName = dbstack('-completenames');
    funcPath = fileparts(funcName(2).file);
    funcName = funcName(2).name;
    
    % Go to the directory containing the file
    currDir = cd(funcPath);
    
    if nargin > 1 && isequal(varargin{1}, 'compile')
        % Standard compilation behaviour
        [varargout{1:nargout}] = ojw_mexcompile(funcName, sourceList, varargin{2:end});
    else
        % Function called without first being compiled
        fprintf('Missing mex file: %s.%s. Will attempt to compile and run.\n', funcName, mexext);
        retval = ojw_mexcompile(funcName, sourceList);
        if retval > 0
            S = which(funcName, '-all'); % Make sure MATLAB registers the new function
            [varargout{1:nargout}] = fevals(funcName, varargin{:});
        else
            % Return to the original directory
            cd(currDir);
            
            % Flag the error
            error('Unable to compile %s.', funcName);
        end
    end
    
    % Return to the original directory
    cd(currDir);
    return
end

% Check for compile flags
for a = nargin:-1:1
    M(a) = varargin{a}(1) == '-';
end

for a = find(~M(:))'
    % Delete current mex
    S = which(varargin{a}, '-all');
    if ~iscell(S)
        S = {S};
    end
    s = cellfun(@(s) strcmp(s(end-numel(varargin{a})-1:end), [varargin{a} '.m']), S);
    if ~any(s)
        error('Function %s not found on the path', varargin{a});
    end
    s = [S{s}(1:end-2) '.' mexext];
    if exist(s, 'file')
        delete(s);
        assert(exist(s, 'file') == 0, 'Could not delete the mex file:\n   %s\nEither it is write protected or it is locked in use by MATLAB.', s);
        S = which(varargin{a}, '-all'); % Refresh the function pointed to
    end
    
    % Compile
    fevals(varargin{a}, 'compile', 0, varargin{M});
    
    % Clear functions and make sure the mex is pointed to
    clear(varargin{a});
    S = which(varargin{a}, '-all');
    assert(exist(s, 'file') ~= 0, 'Compilation of %s failed', varargin{a});
end 
end

%OJW_MEXCOMPILE  Mex compile helper function
%
%   okay = ojw_mexcompile(funcName)
%   okay = ojw_mexcompile(..., inputStr)
%   okay = ojw_mexcompile(..., lastCompiled)
%
% Compile mexed function, given an optional list of source files and other
% compile options. Can optionally check if source files have been modified
% since the last compilation, and only compile if they have.
%
%IN:
%   funcName - string containg the name of the function to compile
%   inputStr - cell array of input strings to be passed to mex, e.g. source
%              file names, include paths, library paths, optimizer flags,
%              etc. Default: {[funcName '.c']}.
%   lastCompiled - datenum of the current mex file. Default: 0 (i.e. force
%                  compilation).
%
%OUT:
%   okay - 1: function compiled; 0: up-to-date, no need to compile; -1:
%          compilation failed.

function okay = ojw_mexcompile(funcName, varargin)

% Determine architecture
is64bit = mexext;
is64bit = strcmp(is64bit(end-1:end), '64');

% Set defaults for optional inputs
sourceList = [funcName '.c'];
lastCompiled = 0;
% Parse inputs
extraOptions = {};
for a = 1:numel(varargin)
    if iscell(varargin{a}) && ischar(varargin{a}{1})
        sourceList = varargin{a};
    elseif isnumeric(varargin{a}) && isscalar(varargin{a})
        lastCompiled = varargin{a};
    elseif ischar(varargin{a})
        extraOptions = [extraOptions varargin(a)];
    end
end
sourceList = [sourceList extraOptions];

okay = 0;
if lastCompiled
    compile = false;
    % Compile if current mex file is older than any of the source files
    % Note: this doesn't consider included files (e.g. header files)
    for a = 1:numel(sourceList)
        dirSource = dir(sourceList{a});
        if ~isempty(dirSource)
            if datenum(dirSource.date) > lastCompiled
                compile = true;
                break;
            end
        end
    end
else
    compile = true;
end

% Exit if not compiling
if ~compile
    return;
end
okay = -1;

debug = false;
cudaOptions = cell(0, 1);
compiler_options = '';
%L = {''};
% Parse the compile options
for a = 1:numel(sourceList)
    if sourceList{a}(1) == '-' % Found an option (not a source file)
        switch sourceList{a}(2)
            case 'N' % Cuda nvcc option
                cudaOptions{end+1} = sourceList{a}(3:end);
                sourceList{a} = '';
            case 'X' % Special option
                [sourceList{a}, co] = feval(sourceList{a}(3:end), debug);
                compiler_options = [compiler_options ' ' co];
            case 'C' % Compiler option
                compiler_options = [compiler_options ' -' sourceList{a}(3:end)];
                sourceList{a} = '';
            case 'g' % Debugging on
                debug = debug | (numel(sourceList{a}) == 2);
        end
    else
        sourceList{a} = ['"' sourceList{a} '"'];
    end
end

L = zeros(numel(sourceList), 1);
gpucc = [];
options_file = [];
% Convert any CUDA files to C++, and any Fortran files to object files
for a = 1:numel(sourceList)
    if isempty(sourceList{a}) || sourceList{a}(1) == '-' % Found an option (not a source file)
        continue;
    end
    [ext, ext, ext] = fileparts(sourceList{a}(2:end-1));
    switch ext
        case '.cu'
            % GPU programming - Convert any *.cu files to *.cpp
            if isempty(gpucc)
                % Create nvcc call
                cudaDir = cuda_path();
                options_file = ['"' tempname '.txt"'];
                fid = fopen(options_file(2:end-1), 'wt'); 
                gpucc = sprintf('"%s%s" --options-file %s', cudaDir, nvcc(), options_file);
                fprintf(fid, ' -D_MATLAB_ -I"%s%sextern%sinclude" -m%d', matlabroot, filesep, filesep, 32*(1+is64bit));
                % Add cuda specific options
                if ~isempty(cudaOptions)
                    fprintf(fid, ' %s', cudaOptions{:});
                end
                if ispc && is64bit
                    cc = mex.getCompilerConfigurations();
                    fprintf(fid, ' -ccbin "%s\\VC\\bin" -I"%s\\VC\\include"', cc.Location, cc.Location);
                end
                % Add any include directories from source list and cuda
                % specific options
                for b = 1:numel(sourceList)
                    if strcmp(sourceList{b}(1:min(2, end)), '-I')
                        fprintf(fid, ' %s', sourceList{b});
                    end
                end
                % Add the debug flag
                if debug
                    fprintf(fid, ' -g -UNDEBUG -DDEBUG');
                else
                    % Apply optimizations
                    fprintf(fid, ' -O3 --use_fast_math');
                end
                fclose(fid);
            end
            % Compile to object file
            outName = ['"' tempname '.o"'];
            cmd = sprintf('%s --compile "%s" --output-file %s', gpucc, sourceList{a}, outName);
            disp(cmd);
            if system(cmd)
                % Quit
                fprintf('ERROR while converting %s to %s.\n', sourceList{a}, outName);
                clean_up([reshape(sourceList(L == 1), [], 1); {options_file}]);
                return;
            end
            sourceList{a} = outName;
            L(a) = 1;
        case {'.f', '.f90'}
            L(a) = 2;
    end
end
% Delete the options file
if ~isempty(options_file)
    delete(options_file(2:end-1));
    options_file = [];
end

% Set the compiler flags
if debug
    flags = '-UNDEBUG -DDEBUG';
else
    flags = '-O -DNDEBUG';
end
if any(L == 1)
    flags = [flags ' ' cuda(debug)];
end
switch mexext
    case {'mexglx', 'mexa64', 'mexmaci64'}
        if ~debug
            compiler_options = [compiler_options ' -O3 -ffast-math -funroll-loops'];
        end
        flags = sprintf('%s CXXOPTIMFLAGS="%s" LDCXXOPTIMFLAGS="%s" LDOPTIMFLAGS="%s"', flags, compiler_options, compiler_options, compiler_options);
    case {'mexw32', 'mexw64'}
        flags = sprintf('%s COMPFLAGS="%s $COMPFLAGS"', flags, compiler_options);
end

% Call mex to compile the code
cmd = sprintf('mex -D_MATLAB_=%d %s%s -output "%s"', [100 1] * sscanf(version(), '%d.%d', 2), flags, sprintf(' %s', sourceList{:}), funcName);
disp(cmd);
try
    eval(cmd);
    okay = 1;
catch me
    fprintf('ERROR while compiling %s\n', funcName);
    fprintf('%s', getReport(me, 'basic'));
end

% Clean up
clean_up(sourceList(L == 1));
end

function clean_up(sourceList)
% Delete the intermediate files
for a = 1:numel(sourceList)
    delete(sourceList{a}(2:end-1));
end
end

function cuda_path_str = cuda_path()
cuda_path_str = user_string('cuda_path');
if ~check_path()
    % Check the environment variables
    cuda_path_str = fullfile(getenv('CUDA_PATH'), '/');
    if check_path()
        user_string('cuda_path', cuda_path_str);
        return;
    end
    % Ask the user to enter the path
    while 1
        cuda_path_str = ask_user_for_directory('CUDA');
        if check_path()
            user_string('cuda_path', cuda_path_str);
            return;
        end
    end
    error('Cuda not found.');
end
% Nested function
    function good = check_path
        % Check the path is valid
        [good, message] = system(sprintf('"%s%s" -h', cuda_path_str, nvcc()));
        good = good == 0;
    end
end

function path_ = nvcc()
path_ = ['bin' filesep 'nvcc'];
if ispc
    path_ = [path_ '.exe'];
end
end

function cuda_sdk_path_str = cuda_sdk
cuda_sdk_path_str = user_string('cuda_sdk_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        cuda_sdk_path_str = [ask_user_for_directory('CUDA SDK') 'C' filesep 'common' filesep];
        if check_path()
             user_string('cuda_sdk_path', cuda_sdk_path_str);
            return;
        end
    end
    warning('Cuda SDK not found.');
end
% Nested function
    function good = check_path
        % Check the path is valid
        good = exist([cuda_sdk_path_str 'inc' filesep 'cutil.h'], 'file');
    end
end

function path_str = ask_user_for_directory(name)
path_str = sprintf('Please select your %s installation directory.', name);
if strncmp(computer, 'MAC', 3) % Is a Mac
    % Give separate warning as the uigetdir dialogue box doesn't have a
    % title
    uiwait(warndlg(path_str))
end
path_str = uigetdir('/', path_str);
if isequal(path_str, 0)
    % User hit cancel or closed window
    error('%s not found.', name);
end
path_str = [path_str filesep];
end

% FEVAL for function which is not a subfunction in this file
function varargout = fevals(func_name_str, varargin)
[varargout{1:nargout}] = feval(str2func(sprintf('@(x) %s(x{:})', func_name_str)), varargin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPECIAL OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the compiler options for cuda
function [str, co] = cuda(debug)
co = '';
% Set gpu code compiler & linker flags
is64bit = mexext;
is64bit = strcmp(is64bit(end-1:end), '64');
bitStr = {'32', '', '64', 'Win32', 'x64'};
cudaDir = cuda_path();
str = sprintf('-I"%sinclude" -L"%slib/%s" -lcudart', cudaDir, cudaDir, bitStr{4+is64bit});
end

function [str, co] = cudasdk(debug)
co = '';
str = sprintf('-I"%sinc"', cuda_sdk);
end

% Create the compiler options for lapack/blas
function [str, co] = lapack(debug)
str = '-lmwlapack -lmwblas'; % Use MATLAB's versions
end

% Create the compiler options for OpenMP
function [str, co] = openmp(debug)
str = '';
if debug
    co = '';
else
    if ispc()
        co = '/openmp';
    else
        co = ''; %'-fopenmp';
    end
end
end

% Create the compiler options for OpenCV
function [str, co] = opencv(debug)
co = '';
opencv_path_str = user_string('opencv_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        path_str = ask_user_for_directory('OpenCV');
        opencv_path_str = [path_str filesep];
        if check_path()
            user_string('opencv_path', opencv_path_str);
            break;
        end
    end
end
str = sprintf('-I"%sinclude/opencv" -L"%slib" -lcv210 -lcvaux210 -lcxcore210', opencv_path_str, opencv_path_str);
% Nested function
    function good = check_path
        % Check the path is valid
        if ispc
            good = exist([opencv_path_str 'cvconfig.h.cmake'], 'file');
        else
            good = exist([opencv_path_str 'cvconfig.h.in'], 'file');
        end
    end
end

% Add the boost library directory
function [str, co] = boost(debug)
co = '';
boost_path_str = user_string('boost_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        boost_path_str = ask_user_for_directory('Boost');
        if check_path()
            user_string('boost_path', boost_path_str);
            break;
        end
    end
end
str = sprintf('-I"%s" -L"%sstage%slib%s"', boost_path_str, boost_path_str, filesep, filesep);
% Nested function
    function good = check_path
        % Check the path is valid
        good = exist(sprintf('%sboost%sshared_ptr.hpp', boost_path_str, filesep), 'file');
    end
end

% Add the directX library directory
function [str, co] = directx(debug)
co = '';
directx_path_str = user_string('directx_sdk_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        directx_path_str = ask_user_for_directory('DirectX SDK');
        if check_path()
            user_string('directx_sdk_path', directx_path_str);
            break;
        end
    end
end
str = sprintf('-L"%sLib%sx%d"', directx_path_str, filesep, 86-22*is64bit());
% Nested function
    function good = check_path
        % Check the path is valid
        if ispc
            good = exist(sprintf('%sLib%sx86%sdxguid.lib', directx_path_str, filesep, filesep), 'file');
        else
            error('DirectX only supported on Windows');
        end
    end
end

% Add the Eigen include directory
function [str, co] = eigen(debug)
co = '';
eigen_path_str = user_string('eigen_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        eigen_path_str = ask_user_for_directory('Eigen');
        if check_path()
            user_string('eigen_path', eigen_path_str);
            break;
        end
    end
end
str = sprintf('-I"%s"', eigen_path_str(1:end-1));
if ~debug
    str = [str ' -DEIGEN_NO_DEBUG'];
end
% Nested function
    function good = check_path
        % Check the path is valid
        good = exist(sprintf('%sEigen%sCore', eigen_path_str, filesep), 'file');
    end
end

% Add the Ceres library directory
function [str, co] = ceres(debug)
co = '';
ceres_path_str = user_string('ceres_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        ceres_path_str = ask_user_for_directory('Ceres');
        if check_path()
            user_string('ceres_path', ceres_path_str);
            break;
        end
    end
end
str = sprintf('-I"%sinclude" -L"%slib" -lceres', ceres_path_str, ceres_path_str);
% Nested function
    function good = check_path
        % Check the path is valid
        good = exist(sprintf('%sinclude%sceres%sceres.h', ceres_path_str, filesep, filesep), 'file');
    end
end

% Add the glog library directory
function [str, co] = glog(debug)
co = '';
glog_path_str = user_string('glog_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        glog_path_str = ask_user_for_directory('glog');
        if check_path()
            user_string('glog_path', glog_path_str);
            break;
        end
    end
end
str = sprintf('-I"%sinclude" -L"%slib" -lglog', glog_path_str, glog_path_str);
% Nested function
    function good = check_path
        % Check the path is valid
        good = exist(sprintf('%sinclude%sglog%slogging.h', glog_path_str, filesep, filesep), 'file');
    end
end

% Add the gflags library directory
function [str, co] = gflags(debug)
co = '';
gflags_path_str = user_string('gflags_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        gflags_path_str = ask_user_for_directory('gflags');
        if check_path()
            user_string('gflags_path', gflags_path_str);
            break;
        end
    end
end
str = sprintf('-I"%sinclude" -L"%slib" -lgflags', gflags_path_str, gflags_path_str);
% Nested function
    function good = check_path
        % Check the path is valid
        good = exist(sprintf('%sinclude%sgflags%sgflags.h', gflags_path_str, filesep, filesep), 'file');
    end
end

% Add the Sophus include directory
function [str, co] = sophus(debug)
co = '';
sophus_path_str = user_string('sophus_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        sophus_path_str = ask_user_for_directory('Sophus');
        if check_path()
            user_string('sophus_path', sophus_path_str);
            break;
        end
    end
end
str = sprintf('-I"%s"', sophus_path_str(1:end-1));
% Nested function
    function good = check_path
        % Check the path is valid
        good = exist(sprintf('%ssophus%sse3.hpp', sophus_path_str, filesep), 'file');
    end
end

% Add the liegroups library
function [str, co] = liegroups(debug)
co = '';
liegroups_path_str = user_string('liegroups_path');
if ~check_path()
    % Ask the user to enter the path
    while 1
        liegroups_path_str = ask_user_for_directory('liegroups');
        if check_path()
            user_string('liegroups_path', liegroups_path_str);
            break;
        end
    end
end
if debug
    str = 'Debug';
else
    str = 'Release';
end
str = sprintf('-I"%s/.." -L"%s/build/%s" -lliegroups', liegroups_path_str, liegroups_path_str, str);
% Nested function
    function good = check_path
        % Check the path is valid
        good = exist(sprintf('%s/se3.hpp', liegroups_path_str), 'file');
    end
end

% Add the ojwul include directory
function [str, co] = ojwul(debug)
co = '';
str = sprintf('-I"%s"', fileparts(fileparts(mfilename('fullpath'))));
end
