%MEX_INTERFACE MATLAB wrapper to an underlying C++ class
classdef mex_interface < handle
    properties (Access = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
        mexHandle; % Handle to the mex function
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = mex_interface(mexfun, varargin)
            this.mexHandle = mexfun;
            this.objectHandle = this.mexHandle('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            this.mexHandle('delete', this.objectHandle);
        end
        
        %% Disp - get the function name
        function disp(this, var_name)
            fprintf('%s is an object instance of %s\n', var_name, func2str(this.mexHandle));
        end

        %% All other methods
        function varargout = subsref(this, s)
            if numel(s) < 2 || ~isequal(s(1).type, '.') || ~isequal(s(2).type, '()')
                error('Not a valid indexing expression')
            end
            [varargout{1:nargout}] = this.mexHandle(s(1).subs, this.objectHandle, s(2).subs{:});
        end
    end
end