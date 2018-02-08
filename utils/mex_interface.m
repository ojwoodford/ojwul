%MEX_INTERFACE MATLAB wrapper to an underlying C++ class
%
% This interface assumes that the mex function uses the following standard
% interface:
%   Construction -    obj = mexfun('new',         ...)
%   Destruction -           mexfun('delete', obj)
%   Other methods - [...] = mexfun('method', obj, ...)
classdef mex_interface < handle
    properties (Access = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
        mexHandle; % Handle to the mex function
    end
    methods
        %% Constructor - Create a new C++ class instance
        % Inputs:
        %    mexfun - handle to the C++ class interface mex.
        %    varargin - arguments passed to the mex when calling 'new'.
        function this = mex_interface(mexfun, varargin)
            this.mexHandle = mexfun;
            this.objectHandle = this.mexHandle('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            if ~isempty(this.objectHandle)
                this.mexHandle('delete', this.objectHandle);
            end
            this.objectHandle = [];
        end
        
        %% Disp - get the function name
        function disp(this, var_name)
            if nargin > 1
                fprintf('%s is an object instance of %s\n', var_name, func2str(this.mexHandle));
            else
                fprintf('Object instance of %s\n', func2str(this.mexHandle));
            end
        end

        %% All other methods
        function varargout = subsref(this, s)
            if numel(s) < 2 || ~isequal(s(1).type, '.') || ~isequal(s(2).type, '()')
                error('Not a valid indexing expression')
            end
            assert(~isempty(this.objectHandle), 'Object not initialized correctly');
            [varargout{1:nargout}] = this.mexHandle(s(1).subs, this.objectHandle, s(2).subs{:});
        end
    end
end