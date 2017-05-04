%SHIFTED_INTERP2 Class for doing shifted linear interpolation multiple times
classdef shifted_interp2 < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = shifted_interp2(im, oobv)
            if nargin < 2
                oobv = NaN;
            end
            this.objectHandle = shifted_interp2_mex(im, oobv);
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            shifted_interp2_mex(this.objectHandle);
        end

        %% Interpolate - Sample the shifted image
        function out = interp2(this, X, Y, varargin)
            if nargin > 4
                % Pass oobv to interpolator
                out = shifted_interp2_mex(this.objectHandle, X, Y, varargin{2});
            else
                out = shifted_interp2_mex(this.objectHandle, X, Y);
            end
        end
    end
end