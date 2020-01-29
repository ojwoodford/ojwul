function estim = cd_learn_normal(varargin)
%CD_LEARN_NORMAL  Demos contrastive divergence learning
%
%   params = cd_learn_normal(options)
%
% Uses contrastive divergence to learn the parameters of a normal
% distribution that training data is generated from, and displays the
% results on completion. Assumes we don't know the normalization factor of
% the model we're learning.
%
% IN:
%   options - Pairs of arguments, the first being a string containing the
%             name of the option, the second being the value of the option.
%             Possible options are:
%      'data' - The number of training data values to use (Default: 1000),
%               or a vector of training data.
%      'niter' - The number of parameter updates to carry out (Default:
%                3000).
%      'mcmc' - The number of MCMC sampling cycles to use (Default: 1).
%      'step' - The step size factor (Default: 0.01).
%      'perturb' - The perturbation distribution from which to draw
%                  updates to the data points in each MCMC cycle (Default:
%                  @(n,sigma) randn(n,1)*sigma).
%      'start' - The parameters to initialize our proposal distribution to,
%                (Default: [0 1]).
%      'finish' - The parameters of the distribution we are aiming for
%                 (Default: [2*rand(1) 0.5+rand(1)]). Ignored if training
%                 data is given.
%      'show_results' - Results graphically displayed if non-zero (Default:
%                       true).
%
% OUT:
%   params - (niter+1)x2 matrix of model parameters at each iteration.

% Load options
opts.data = 1000;
opts.niter = 3000;
opts.mcmc = 1;
opts.step = 0.01;
opts.start = [0 1];
opts.finish = [2*rand(1) 0.5+rand(1)];
opts.perturb = @(n,sigma) randn(n, 1) * sigma;
opts.show_output = true;
opts = argparse(opts, varargin{:});

% Generate the data and our initial state
if numel(opts.data) == 1
    X = randn(opts.data, 1);
    X = (X * opts.finish(2)) + opts.finish(1);
else
    X = opts.data(:);
    opts.data = numel(opts.data);
end
params = opts.start(:);

% Buffers for recording display data
estim = zeros(2, opts.niter+1);
estim(:,1) = params;
compare = estim;

for iter = 1:opts.niter    
    % Create confabulations using MCMC (Metropolis algorithm)
    x = X;
    px = calc_unnormalized_probability(x, params);
    for k = 1:opts.mcmc
        % Sample from our perturbation distribution
        y = x + opts.perturb(opts.data, params(2));
        
        % Reject some of the samples
        py = calc_unnormalized_probability(y, params);
        L = py > (px .* rand(opts.data, 1));
        x(L) = y(L);
        px(L) = py(L);
    end

    % Calculate the gradient of our energy function
    grad = opts.step * (calc_expectation_logfx(X, params) - calc_expectation_logfx(x, params));

    % Update our parameters
    params = params + grad;
    if params(2) <= 0
        params(2) = 1e-7; % Ensure sd is positive
    end
    
    % Record this for display later
    estim(:,iter+1) = params;
end

% Display the results
estim = estim(:,1:iter+1)';
if opts.show_output
    mu = mean(X);
    sd = sqrt(var(X, 1));
    clf
    subplot(311);
    plot([0 iter], [mu mu], 'r-');
    hold on
    plot(0:iter, estim(:,1), 'b-');
    data.l1 = axis;
    data.l1 = data.l1(3:4);
    data.h1 = plot([0 0], data.l1, 'k-');
    title(sprintf('Mean - actual: %g, learnt: %g', mu, estim(end,1)));
    legend('Linear solve', 'Contrastive divergence', 'Location', 'Best');
    subplot(312);
    plot([0 iter], [sd sd], 'r-');
    hold on
    plot(0:iter, estim(:,2), 'b-');
    data.l2 = axis;
    data.l2 = data.l2(3:4);
    data.h2 = plot([0 0], data.l2, 'k-');
    title(sprintf('Standard deviation - actual: %g, learnt: %g', sd, estim(end,2)));
    legend('Linear solve', 'Contrastive divergence', 'Location', 'Best');
    subplot(313);
    plot(X, abs(randn(size(X)))*0.1+0.01, 'k.');
    hold on
    data.normal{1} = [-3:0.25:-0.75 -0.5:0.1:0.5 0.75:0.25:3];
    data.normal{2} = exp(-(data.normal{1} .^ 2));
    [X, Y] = calc_normal_dist(mu, sd, data.normal);
    plot(X, Y, 'r-');
    [X, Y] = calc_normal_dist(estim(1,1), estim(1,2), data.normal);
    data.h3 = plot(X, Y, 'b-');
    X = axis;
    axis([X(1)-1 X(2)+1 X(3) X(4)]);
    title('Estimated distributions - use mouse or arrow keys to select a column above');
    legend('Data points', 'Linear solve', 'Contrastive divergence', 'Location', 'Best');
    % Initialize the callbacks
    data.col = 0;
    data.niter = iter;
    data.step = round(iter/300);
    data.estim = estim;
    data.compare = compare;
    set(gcf, 'WindowButtonDownFcn', @mousedown_callback, ...
        'KeyPressFcn', @keypress_callback, ...
        'UserData', data);
end
if nargout == 0
    clear estim
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad = calc_expectation_logfx(X, params)
A = X - params(1);
B = (A .^ 2) / (params(2) ^ 3);
A = A / (params(2) ^ 2);
grad = [mean(A); mean(B)];
end

function px = calc_unnormalized_probability(X, params)
px = exp(-((X - params(1)) .^ 2)  / (2 * (params(2) ^ 2)));
end

% The functions below are only used for displaying results
function [X, Y] = calc_normal_dist(mu, sd, range)
X = range{1} * sd + mu;
Y = range{2} / (sd * sqrt(2 * pi));
end

function data = plot_data(data)
data.col = mod(data.col, data.niter);

% Plot the new lines
set(data.h1, 'XData', data.col([1 1]), 'YData', data.l1);
set(data.h2, 'XData', data.col([1 1]), 'YData', data.l2);
[X, Y] = calc_normal_dist(data.estim(data.col+1,1), data.estim(data.col+1,2), data.normal);
set(data.h3, 'XData', X, 'YData', Y);
drawnow();
end

function keypress_callback(fig, event_data)
data = get(fig, 'UserData');
key = lower(event_data.Character);
% Which keys do what
switch key
    case 28 % Left
        data.col = data.col - data.step;
    case 29 % Right
        data.col = data.col + data.step;
    case 31 % Down
        data.col = data.col - 10 * data.step;
    case 30 % Up
        data.col = data.col + 10 * data.step;
end
data = plot_data(data);
set(fig, 'UserData', data);
end

function mousedown_callback(fig, event_data)
% Detect which function is required
button = get(fig, 'SelectionType');
data = get(fig, 'UserData');
switch button
    case 'extend'
    case 'alt'
    otherwise
        % Button 1:
        data.col = get(gca, 'CurrentPoint');
        data.col = round(data.col(1));
        data = plot_data(data);
        set(fig, 'UserData', data);
end
end

function opts = argparse(opts, varargin)
if isempty(varargin)
    inopts = struct([]);
else
    inopts = struct(varargin{:});
end

fn = fieldnames(inopts);
for i = 1:length(fn)
    if isfield(opts, fn{i})
        opts.(fn{i}) = inopts.(fn{i});
    else
        error('Bad argument: ''%s''', fn{i});
    end
end
end
