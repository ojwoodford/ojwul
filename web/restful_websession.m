%RESTFUL_WEBSESSION Maintains cookies across requests within a web session
% 
% This class provides a RESTful API for accessing web resources, that
% maintains cookies across requests. This can be useful for sessions that
% require authentication. It is built atop MATLAB's RESTful http interface.
%    
%Example:
%   h = restful_websession('https://bbc.com/news');
%   response = h.get();
%
%   See also MATLAB.NET.HTTP.REQUESTMESSAGE.

% Copyright (C) Oliver Woodford 2025

classdef restful_websession < handle
    properties (Hidden = true, SetAccess = protected)
        uri;     % matlab.net.URI
        request; % matlab.net.http.RequestMessage
        options; % matlab.net.http.HTTPOptions
    end

    methods
        % Constructor - takes in the url we'll be communicating with
        function this = restful_websession(url)
            this.uri = matlab.net.URI(url);
            this.request = matlab.net.http.RequestMessage();
            this.options = matlab.net.http.HTTPOptions('ProgressMonitorFcn', @web_progressbar, 'ResponseTimeout', 20, 'DataTimeout', 20);
        end

        % SET_OPTIONS - takes in Name,Value pairs
        function set_options(this, varargin)
            for a = 1:2:numel(varargin)
                this.options.(varargin{a}) = varargin{a+1};
            end
        end

        % GET - takes in query Name,Value pairs
        function response = get(this, varargin)
            response = get_helper(this, {}, varargin);
        end

        % DOWNLOAD - get request that saves to a file
        function response = download(this, fname, varargin)
            % Construct the consumer
            consumer = matlab.net.http.io.FileConsumer(fname, 'w', 'n');
            % Send the request
            response = get_helper(this, {consumer}, varargin{:});
        end
    end

    methods (Hidden = true, Access = protected)
        function response = get_helper(this, consumer, varargin)
            % Construct the query
            this.request.Method = matlab.net.http.RequestMethod.GET;
            this.uri.Query = matlab.net.QueryParameter(varargin{1:end});
            % Send the request
            response = request_set_cookies(this, consumer);
            % Check for a 303 redirect (which matlab.net.http.RequestMessage doesn't handle)
            if response.StatusCode == matlab.net.http.StatusCode.SeeOther
                % Redirect
                this.uri = matlab.net.URI(response.Header(arrayfun(@(s) strcmp(s.Name, 'Location'), response.Header)).Value);
                % Send the new request
                response = request_set_cookies(this, consumer);
            end
        end

        function response = request_set_cookies(this, consumer)
            % Send the request
            [response, ~, history] = this.request.send(this.uri, this.options, consumer{:});
            % Update the URI
            this.uri = history(end).URI;
            % Store cookies for future requests
            for record = history
                cookieFields = record.Response.getFields('Set-Cookie');
                if isempty(cookieFields)
                    continue;
                end
                this.request.replaceFields(matlab.net.http.field.CookieField([cookieFields.convert().Cookie]));
            end
        end
    end
end
