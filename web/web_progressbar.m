%WEB_PROGRESSBAR A progress bar for web downloads
% 
% This class works with MATLAB's RESTful http interface to provide a
% visual progress bar showing the progress of uploads and downloads.
%    
%Example:
%   import matlab.net.http.*
%   opt = HTTPOptions('ProgressMonitorFcn', @web_progressbar, 'UseProgressMonitor', true);
%   response = send(RequestMessage(), 'https://bbc.com/news', opt)
%
%   See also MATLAB.NET.HTTP.PROGRESSMONITOR, OJW_PROGRESSBAR.

% Copyright (C) Oliver Woodford 2025

classdef web_progressbar < matlab.net.http.ProgressMonitor
    properties
        PBH
        Value
        Direction
        oldDirection
        oldValue
        oldMax
    end
    
    methods
        function this = web_progressbar
            this.Interval = 0.5;
        end
        
        function done(this)
            this.close();
        end
        
        function delete(this)
            this.close();
        end
        
        function set.Value(this, value)
            this.Value = value;
            this.update();
        end
    end
    
    methods (Access = private)
        function update(this)
            import matlab.net.http.*
            % Check for changes in status
            if isempty(this.Value) || isempty(this.Max) || isempty(this.Direction)
                close(this);
                this.oldValue = 0;
                return;
            end
            if ~isempty(this.PBH) && (this.Value < this.oldValue || ~isequal(this.Max, this.oldMax) || ~isequal(this.Direction, this.oldDirection))
                close(this);
            end
            this.oldValue = this.Value;
            this.oldMax = this.Max;
            this.oldDirection = this.Direction;
            if this.Max < 1e5
                % Do not display for transfers below 100KB
                return;
            end

            % Create or update the progress bar
            if isempty(this.PBH)
                if this.Direction == MessageType.Request
                    msg = 'Uploading...';
                else
                    msg = 'Downloading...';
                end
                this.PBH = ojw_progressbar(msg, this.Value, this.Max, 0.2);
            else
                update(this.PBH, this.Value);
            end
        end

        function close(this)
            if isempty(this.PBH)
                return;
            end
            close(this.PBH);
            this.PBH = [];
        end
    end
end
