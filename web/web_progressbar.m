classdef web_progressbar < matlab.net.http.ProgressMonitor
    properties
        ProgHandle
        Value uint64
        Direction matlab.net.http.MessageType
        oldDirection matlab.net.http.MessageType
        oldValue uint64
        oldMax uint64
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
                this.oldMax = [];
                this.oldDirection = [];
                return;
            end
            if ~isempty(this.ProgHandle) && (this.Value < this.oldValue || ~isequal(this.Max, this.oldMax) || ~isequal(this.Direction, this.oldDirection))
                close(this);
            end
            this.oldValue = this.Value;
            this.oldMax = this.Max;
            this.oldDirection = this.Direction;
            if this.Max < 1e5
                % Do not display for small transfers
                return;
            end

            % Create or update the progress bar
            if isempty(this.ProgHandle)
                if this.Direction == MessageType.Request
                    msg = 'Uploading...';
                else
                    msg = 'Downloading...';
                end
                this.ProgHandle = ojw_progressbar(msg, this.Value, this.Max, 0.2);
            else
                update(this.ProgHandle, this.Value);
            end
        end

        function close(this)
            if ~isempty(this.ProgHandle)
                close(this.ProgHandle);
                this.ProgHandle = [];
            end
        end
    end
end
