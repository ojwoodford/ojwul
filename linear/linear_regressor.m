classdef linear_regressor < handle
    properties (SetAccess = private, Hidden = true)
        parameters;
        regularization_lambda = 0;
        isquadratic = false;
    end
    methods
        function this = linear_regressor(quad, reg)
            if nargin > 0
                this.isquadratic = quad;
                if nargin > 1
                    this.regularization_lambda = reg;
                end
            end
        end
        
        function train(this, X, y)
            % Preprocess the data
            X = quadraticize(this, X');
            norm_mat = compute_normalization(X);
            X(:,end+1) = 1;
            X = X * norm_mat;
            
            % Train
            if this.regularization_lambda
                % Minimize the regularized cost using a local optimizer
                this.parameters = fmincg(@(theta) linear_regression_cost(X, y, theta, this.regularization_lambda), zeros(size(X, 2), 1), optimset('MaxIter', 200, 'GradObj', 'on'));
            else
                % Solve the normal equation linear system
                this.parameters = (X' * X) \ (X' * y(:));
            end
            
            % Apply the normalization
            this.parameters = norm_mat * this.parameters;
        end
        
        function y = test(this, X)
            % Preprocess the data
            X = quadraticize(this, X');
            X(:,end+1) = 1;
            
            % Compute the output
            y = X * this.parameters;
        end
        
        function X = quadraticize(this, X)
            if this.isquadratic
                [m, n] = size(X);
                X = [X reshape(bsxfun(@times, X, reshape(X, m, 1, n)), m, n*n)];
            end
        end
    end
end

function norm_mat = compute_normalization(X)
m = mean(X, 1);
X = bsxfun(@minus, X, m);
s = 1 ./ (sqrt(mean(X .* X, 1)) + 1e-38);
norm_mat = diag(s);
norm_mat(end+1,:) = -(m .* s);
norm_mat(end,end+1) = 1;
end

function [J, grad] = linear_regression_cost(X, y, theta, lambda)
if nargin < 4
    lambda = 0;
end
m = numel(y); % number of training examples
hx = X * theta;
grad = hx - y(:);
J = grad' * grad;
if lambda
    J = J + (theta(2:end)' * theta(2:end)) * lambda;
end
J = J / (2 * m);
if nargout > 1
    grad = sum(bsxfun(@times, grad, X), 1)' / m;
    if lambda
        grad(2:end) = grad(2:end) + theta(2:end) * (lambda / m);
    end
end
end

    