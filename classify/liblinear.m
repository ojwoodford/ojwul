%LIBLINEAR Interface to the liblinear library through a standard API
%
% Methods:
%   obj = liblinear(); % Constructor - does nothing
%   train(obj, data, labels, options);
%   labels = train(obj, data);
%      data - MxN array of N features of M dimensions.
%      labels - 1xN logical vector of input classes or output class
%               probabilities.
%      options - structure containing the following fields:
%


classdef liblinear < handle
    properties (Access = private, Hidden = true)
        params;
    end
    
    methods
        % Constructor - does nothing
        function this = liblinear()
        end

        % Train
        function train(this, data, labels, varargin)            
            % Train the SVM
            optString = this.construct_options(varargin{:});
            
            % Add weights (only for C-SVC but no harm adding anyway)
            labels = labels(:)*2-1;
            w = hist(labels, [-1 1]);
            w = w ./ sum(w);
            optString = [optString sprintf(' -w-1 %d -w1 %d', w)];
            
            this.params = linear_train(labels, sparse(data), optString, 'col');
        end
        
        % Test
        function [probs, scores, accuracy] = test(this, data)
            % Run the classifier
            [scores, accuracy, probs] = linear_predict(ones(size(data, 2), 1), data', get_params(this));
        end
        
        % Perform cross validation to find the best parameters
        function bestParams = param_search(this, data, labels, varargin)
            % Grid search here depends on type, only c, or c and gamma
            if isfield(varargin{1},'cSearch')
                cSearch = varargin{1}.cSearch;
            else
                cSearch =  logspace(-10,1,10);
            end
            
            % Ensure we have some cross val folds
            if ~isfield(varargin{1},'val') || varargin{1}.val < 2
                varargin{1}.val = 3;
            end
            
            % Perform grid search
            accuracy = zeros(numel(cSearch), 1);
            for c = 1:numel(cSearch)
                cost = cSearch(c);
                varargin{1}.cost = cost;
                train(this, data, labels, varargin{1});
                accuracy(c) = this.params;
            end
            [~, bestCI] = max(accuracy);
            bestParams.cost = cSearch(bestCI);
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
        % options:
        % -s type : set type of solver (default 1)
        %   for multi-class classification
        % 	 0 -- L2-regularized logistic regression (primal)
        % 	 1 -- L2-regularized L2-loss support vector classification (dual)
        % 	 2 -- L2-regularized L2-loss support vector classification (primal)
        % 	 3 -- L2-regularized L1-loss support vector classification (dual)
        % 	 4 -- support vector classification by Crammer and Singer
        % 	 5 -- L1-regularized L2-loss support vector classification
        % 	 6 -- L1-regularized logistic regression
        % 	 7 -- L2-regularized logistic regression (dual)
        %   for regression
        % 	11 -- L2-regularized L2-loss support vector regression (primal)
        % 	12 -- L2-regularized L2-loss support vector regression (dual)
        % 	13 -- L2-regularized L1-loss support vector regression (dual)
        % -c cost : set the parameter C (default 1)
        % -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
        % -e epsilon : set tolerance of termination criterion
        % 	-s 0 and 2
        % 		|f'(w)|_2 <= eps*min(pos,neg)/l*|f'(w0)|_2,
        % 		where f is the primal function and pos/neg are # of
        % 		positive/negative data (default 0.01)
        % 	-s 11
        % 		|f'(w)|_2 <= eps*|f'(w0)|_2 (default 0.001) 
        % 	-s 1, 3, 4 and 7
        % 		Dual maximal violation <= eps; similar to libsvm (default 0.1)
        % 	-s 5 and 6
        % 		|f'(w)|_inf <= eps*min(pos,neg)/l*|f'(w0)|_inf,
        % 		where f is the primal function (default 0.01)
        % 	-s 12 and 13\n"
        % 		|f'(alpha)|_1 <= eps |f'(alpha0)|,
        % 		where f is the dual function (default 0.1)
        % -B bias : if bias >= 0, instance x becomes [x; bias]; if < 0, no bias term added (default -1)
        % -wi weight: weights adjust the parameter C of different classes (see README for details)
        % -v n: n-fold cross validation mode
        % -q : quiet mode (no outputs)
        function str = construct_options(opts)
            % Defaults
            c = 1; % Cost (regulariser)
            s = 1; % Type
            v = 0; % Cross validation folds (0 = no cross val)
            q = 1; % Quiet mode
            if nargin > 0
                % Type
                if isfield(opts, 'type')
                    if ~ismember([0:7 11:13], opts.type)
                        error('Invalid type.');
                    else
                        s = opts.type;
                    end
                end
                % Cost
                if isfield(opts, 'cost')
                    c = opts.cost;
                end
                % Cross validation folders
                if isfield(opts, 'val')
                    v = opts.val;
                end
                % Quiet mode
                if isfield(opts, 'quiet')
                    q = opts.quiet ~= 0;
                end
            end
            str = sprintf('-c %d -s %d', c, s);
            
            if v > 0
                str = [str sprintf(' -v %d', v)];
            end
            
            if q == 1
                str = [str ' -q'];
            end
        end
    end
end