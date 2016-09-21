%LIBSVM Interface to the libsvm library through a standard API
%
% Methods:
%   obj = libsvm(options); % Constructor - does nothing
%   train(obj, data, labels);
%   labels = train(obj, data);
%      data - MxN array of N features of M dimensions.
%      labels - 1xN logical vector of input classes or output class
%               probabilities.
%      options - structure containing the following fields:
%


classdef libsvm < handle
    properties (Access = private, Hidden = true)
        params;
        optString;
    end
    
    methods
        % Constructor - set the SVM parameters
        function this = libsvm(varargin)
            this.optString = this.construct_options(varargin{:});
        end

        % Train
        function train(this, data, labels)
            % Add weights (only for C-SVC but no harm adding anyway)
            labels = labels(:)*2-1;
            w = hist(labels, [-1 1]);
            w = w ./ sum(w);
            
            % Train the SVM
            this.params = svm_train(labels, data', [this.optString sprintf(' -w-1 %d -w1 %d', w)]);
        end
        
        % Test
        function [probs, scores, accuracy] = test(this, data)
            % Run the classifier
            [scores, accuracy, probs] = svm_predict(ones(size(data, 2), 1), data', get_params(this));
        end
        
        % Perform cross validation to find the best parameters
        function bestParams = param_search(this, data, labels, varargin)
            % Grid search here depends on type, only c, or c and gamma
            warning(0, 'Param_search unimplemented! Defaults will be used unless specified.');
            bestParams.cost = 1;
        end
        
        % Parameter getting and setting, for precomputing models
        function params = get_params(this)
            if isempty(this.params)
                error('SVM has not been trained!');
            end
            params = this.params;
        end
        function set_params(this, params)
            this.params = params;
        end
        function save_params(this, fname)
            svm_params = get_params(this);
            save(fname, 'svm_params');
        end
        function load_params(this, fname)
            set_params(this, getfield(load(fname), 'svm_params'));
        end
    end
    
    methods (Static = true, Hidden = true, Access = private)
        % Turn the standard options structure into a string
        function str = construct_options(opts)
            % Defaults
            s = 0; % C-SVC type
            t = 0; % Linear kernel
            b = 0; % Not probabilistic
            c = 1; % Cost (regulariser)
            g = []; % Gamma in kernel func
            d = 3; % Degree in kernel func
            r = 0; % coef0 in kernel func
            n = 0.5; % Nu parameter
            if nargin > 0
                % SVM type
                if isfield(opts, 'type')
                    s = find(ismember({'c-svc', 'nu-svc', 'one-class', 'epsilon-svr', 'nu-svr'}, lower(opts.type)), 1) - 1;
                    if isempty(s)
                        error('Type name not recognized.');
                    end
                end
                % Kernel type
                if isfield(opts, 'kernel')
                    t = find(ismember({'linear', 'polynomial', 'radial', 'sigmoid', 'intersection', 'precomp'}, opts.kernel), 1) - 1;
                    if isempty(t)
                        error('Kernel name not recognized.');
                    end
                end
                % Probabilistic
                if isfield(opts, 'probabilistic')
                    b = opts.probabilistic ~= 0;
                end
                % Cost
                if isfield(opts, 'cost')
                    c = opts.cost;
                end
                % Gamma
                if isfield(opts, 'gamma')
                    g = opts.gamma;
                end
                % Degree
                if isfield(opts, 'degree')
                    d = opts.degree;
                end
                % coef0
                if isfield(opts, 'coef0')
                    r = opts.coef0;
                end
                % Nu
                if isfield(opts, 'nu')
                    n = opts.coef0;
                end
            end
            if ~isempty(g)
                str = sprintf('-c %d -g %d -s %d -t %d -b %d -d %d -r %d -n %d', c, g, s, t, b, d, r, n);
            else
                str = sprintf('-c %d -s %d -t %d -b %d -d %d -r %d -n %d', c, s, t, b, d, r, n);
            end
        end
    end
end