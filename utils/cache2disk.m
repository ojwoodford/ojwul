%CACHE2DISK Wrapper class for caching slow-to-load data on disk
%
%    buffer = cache2disk(cache_file, load_func)
%    obj = buffer{key}
%
% This class implements a cache on disk, which can improve efficiency when
% loading slow-to-load objects several times, across different MATLAB
% sessions.
%
% IN:
%   cache_file - Filename or full path to the cache .mat file.
%   load_func - Handle to a function which takes a scalar key as input and
%               returns an object.
%   key - scalar key which is passed to read_fun to load an object.
%
% OUT:
%   buffer - handle to the cache.
%   obj - cached object.

% Copyright (C) Oliver Woodford 2015

classdef cache2disk < handle
    properties (Hidden = true, SetAccess = private)
        load_func; % Function to read in an object
        cache_handle; % The matfile object to use as the cache
    end
    
    methods
        % Constructor
        function this = cache2disk(cache_file, load_fun)
            this.load_func = load_fun;
            this.cache_handle = matfile(cache_file, 'Writable', true);
        end
        % The main function - get
        function A = get(this, n)
            if nargin < 2 || ~isscalar(n)
                error('Only one object can be got at a time');
            end
            key = sprintf('obj_%.16g', n);
            try
                A = this.cache_handle.(key);
            catch
                A = this.load_func(n);
                this.cache_handle.(key) = A;
            end
        end
        % Forward calls like cache(a) to get
        function A = subsref(this, frame)
            switch frame(1).type
                case {'()', '{}'}
                    if numel(frame(1).subs) ~= 1
                        error('Only one dimensional indexing supported');
                    end
                    A = get(this, frame(1).subs{1});
                case '.'
                    if strcmp(frame(1).subs, 'get')
                        % Forward these references to the relevant method
                        A = builtin('subsref', this, frame);
                    else
                        error('%s is not a public property or method of the cache class.', frame(1).subs);
                    end
            end
        end
    end
end
