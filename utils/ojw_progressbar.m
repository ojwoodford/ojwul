%OJW_PROGRESSBAR  Simple progress bar implementation
%
%   [this, retval] = ojw_progressbar(tag, proportion[, min_update_interval])
%
% Starts, updates and closes a progress bar according to the proportion of
% time left.
%
% IN:
%   tag - String that appears on progress bar, specific to each function
%         calling OJW_PROGRESSBAR.
%   proportion - Proportion of time elapsed, between 0 and 1.
%   min_update_interval - Minimum time (in seconds) between updates of the
%                         progress bar. The value is kept while the bar is
%                         alive. Default: 0.5.
%
% OUT:
%   this - Handle to the progressbar object.
%   retval - 2 iff time since last report > min_update_interval,
%            1 iff progress bar initialized or reset,
%            0 otherwise.

% Based on Andrew Fitzgibbon's awf_progressbar

classdef ojw_progressbar < handle
    properties (Hidden = true, SetAccess = private)
        bar;
        min_update;
        prop;
        start_prop;
        next_prop;
        timer;
        last_update;
        text_version;
        tag;
        tag_title;
        total;
        inverse_total;
    end
    methods
        function [this, retval] = ojw_progressbar(tag, proportion, total, min_update_interval)
            % Check the input arguments
            if nargin < 2
                error('At least 2 input arguments expected');
            end
            if ~ischar(tag)
                error('First argument should be a string');
            end
            if ~isscalar(proportion) || proportion < 0
                error('Second argument should be a non-negative scalar');
            end
            if nargin > 2
                if ~isscalar(total) || total < 0
                    error('Third argument should be a non-negative scalar');
                end
                total = double(total);
            else
                total = 1; % Default total proportion
            end
            if nargin > 3
                if ~isscalar(min_update_interval) || min_update_interval < 0
                    error('Fourth argument should be a non-negative scalar');
                end
                min_update_interval = double(min_update_interval);
            else
                min_update_interval = 0.5; % Default seconds between updates
            end

            % Record the time
            curr_time = clock();
            proportion = double(proportion); % Must be a double

            % Check the tag exists
            tag_title = tag;
            tag = tag(isstrprop(tag, 'alphanum'));
            v = ojw_progressbar.GetSetPersistent(tag);
            if isempty(v)
                % No tag by this name
                if proportion >= total
                    % No need to create one
                    retval = 0;
                    return;
                end

                % Create a data structure for this tag
                this.bar = [];
                this.min_update = min_update_interval; % Default seconds between updates
                this.prop = proportion;
                this.start_prop = proportion;
                this.next_prop = proportion;
                this.timer = curr_time;
                this.last_update = curr_time;
                this.tag_title = tag_title;
                this.tag = tag;
                this.total = total;
                this.inverse_total = 1 / total;
                this.text_version = ojw_progressbar.GetSetPersistent('text_version');
                
                % Update our global variable with the changes to this tag
                ojw_progressbar.GetSetPersistent(tag, this);
                
                % Ensure it gets closed when the function exits
                [~, varname] = fileparts(tempname());
                assignin('caller', varname, onCleanup(@() ojw_progressbar(tag, 1)));
            else
                % Use the cached the data structure
                this = v;
            
                % Update the total if a new one is given
                if nargin > 2
                    this.total = total;
                    this.inverse_total = 1 / total;
                    this.next_prop = this.prop;
                    % Update the minimum update interval if a new one is given
                    if nargin > 3
                        this.min_update = min_update_interval;
                    end
                end
            end

            % Update the bar
            retval = update(this, proportion, curr_time); 
        end
        
        % Function for direct updating
        function retval = update(this, proportion, curr_time)
            % Check for a fast (i.e. no) update
            if proportion < this.next_prop && proportion >= this.prop
                retval = 1;
                return;
            end

            % Check for an end to the progress bar
            retval = 0;
            if proportion >= this.total
                % Close the progress bar
                if this.text_version
                    fprintf([repmat(' ', 1, 200) repmat('\b', 1, 200)]);
                else
                    close(this.bar);
                    drawnow();
                end
                ojw_progressbar.GetSetPersistent(this.tag, []);
                this.prop = this.total;
                return;
            end
            
            if nargin < 3
                curr_time = clock();
                proportion = double(proportion); % Must be a double
            end

            % Check to see if we haven't started again
            if proportion < this.prop || proportion == this.start_prop
                retval = 1;
            elseif etime(curr_time, this.last_update) >= this.min_update
                % An update of the progress bar is required   
                retval = 1 + (proportion > this.start_prop);
            else
                % No update required
                return;
            end

            this.prop = proportion;
            proportion = proportion * this.inverse_total;
            if retval == 2
                this.last_update = curr_time;
                t_elapsed = etime(curr_time, this.timer);
                prop_done = this.prop - this.start_prop;
                t_remaining = ((this.total - this.prop) * t_elapsed) / prop_done;
                this.next_prop = min(this.prop + prop_done * this.min_update * 1.1 / t_elapsed, this.total);
                newtitle = sprintf('Elapsed: %s', timestr(t_elapsed));
                if proportion > 0.01 || t_elapsed > 30
                    if t_remaining < 600
                        newtitle = sprintf('%s, Remaining: %s          ', newtitle, timestr(t_remaining));
                    else
                        newtitle = sprintf('%s, ETA: %s', newtitle, datestr(datenum(curr_time) + (t_remaining * 1.15741e-5), 0));
                    end
                end
            else
                % Reset the information
                this.start_prop = this.prop;
                this.next_prop = this.prop;
                this.timer = curr_time;
                this.last_update = this.timer;
                newtitle = 'Starting...';
            end

            % Update the waitbar
            if this.text_version
                % Text version
                proportion = floor(proportion * 50);
                str = sprintf('  %s |%s%s| %s', this.tag_title, repmat('#', 1, proportion), repmat(' ', 1, 50 - proportion), newtitle);
                fprintf([str repmat('\b', 1, numel(str))]);
            else
                % Graphics version
                if ishandle(this.bar)
                    waitbar(proportion, this.bar, newtitle);
                else
                    % Create the waitbar
                    this.bar = waitbar(proportion, newtitle, 'Name', this.tag_title);
                end
            end
            drawnow();
        end
    end
    methods (Static = true, Access = private)
        function val = GetSetPersistent(tag, val)
            % Ensure the global data structure exists
            persistent ojw_progressbar_data
            if isempty(ojw_progressbar_data)
                ojw_progressbar_data.text_version = ~usejava('desktop');
            end
            if nargin < 2
                if isfield(ojw_progressbar_data, tag)
                    val = ojw_progressbar_data.(tag);
                else
                    val = [];
                end
            else
                ojw_progressbar_data.(tag) = val;
            end
        end
    end
end

% Time string function
function str = timestr(t)
s = rem(t, 60);
m = rem(floor(t/60), 60);
h = floor(t/3600);

if h > 0
    str= sprintf('%dh%02dm%02.0fs', h, m, s);
elseif m > 0
    str = sprintf('%dm%02.0fs', m, s);
else
    str = sprintf('%2.1fs', s);
end
end
