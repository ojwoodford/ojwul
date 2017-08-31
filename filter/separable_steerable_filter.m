%SEPARABLE_STEERABLE_FILTER Compute the basis and weights of separable steerable filters
%
%   [basis, weight_fun] = separable_steerable_filter(r, coeff, X)
%
% Decompose a steerable filter formed from an odd or even parity polynomial
% times a Gaussian window into a separable basis using the equations in
% Appendix D of:
% "The Design and Use of Steerable Filters", Freeman & Adelson, TPAMI 1991
%
%IN:
%   r - Scalar radius (s.d.) of the Gaussian window.
%   coeff - 1xN vector of coefficients of the polynomial, highest order first.
%   X - 1xM vector of x (and y) locations to sample the filter at.
%
%OUT:
%   basis - Mx2xN array of N pairs of separable filters (y direction filter first).
%   weight_fun - handle to a function which takes an angle as input and
%                returns the N weights for the basis filter responses which
%                should be summed to give the output of the original filter
%                at that angle.

function [basis, weight_fun] = separable_steerable_filter(r, coeff, X)

% Default parameters for visualizing
if nargin < 3
    if nargin < 2
        if nargin < 1
            % Radius of the Gaussian mask
            r = 1;
        end
        % Coefficients of the polynomial, highest order first
        coeff = [4 0 -2];
    end
    % Values of x (and y) to evaluate the filter basis at
    X = linspace(-r*4, r*4, 400);
end
X = reshape(X, 1, []);

% Check filter is odd or even parity
assert(all(coeff(1:2:end) == 0) || all(coeff(2:2:end) == 0), 'Filter not odd or even parity');

% Construct the filter function
f = @(theta, x, y) filter_val(r, coeff, theta, x, y);

% Compute regularly spaced angles
N = find(fliplr(coeff), 1, 'last') - 1;
theta = linspace(0, pi, N+2);
theta = theta(1:end-1);

% Compute the filters at those angles
F = f(reshape(theta, 1, 1, 1, []), X, X');

% Compute the weight function of eqn. 42
K = (0:N)';
weight_fun = @(theta) bsxfun(@times, (-1 .^ K) .* (factorial(N) ./ (factorial(K) .* factorial(N - K))), bsxfun(@power, cos(theta), N-K) .* bsxfun(@power, sin(theta), K));

% Compute the K matrix of eqn. 43
K = weight_fun(theta)';

% Compute the basis filters
basis = K \ reshape(permute(F, [4 1 2 3]), numel(theta), []);
M = numel(X);
basis = reshape(basis', M, M, N+1);

% Compute the separable x and y directions of each filter
[x, x] = max(abs(reshape(basis, M*M, N+1)), [], 1);
[y, x] = ind2sub([M M], x);
B = basis;
basis = zeros(M, 2, N+1);
for a = 1:N+1
    v = 1 ./ sqrt(abs(B(y(a),x(a),a)));
    basis(:,1,a) = B(:,x(a),a) * v * sign(B(y(a),x(a),a));
    basis(:,2,a) = B(y(a),:,a) * v;
end

% Display the filters and basis
if nargout == 0
    figure(1); clf reset; sc(F, 'diff');
    figure(2); clf reset; sc(reshape(B, M, M, 1, N+1), 'diff');
    clear basis weight_fun
end
end

% Function to construct the Gaussian mask
function g = gauss_mask(x, y, r)
g = exp(-0.5 * bsxfun(@plus, x .* x, y .* y) / (r * r)) / (r * sqrt(2 * pi));
end

% Function to construct the oriented filter
function v = filter_val(r, coeff, theta, x, y)
% Compute x' of eqn. 38
v = bsxfun(@plus, bsxfun(@times, x, cos(theta)), bsxfun(@times, y, sin(theta)));
% Multiply the Gaussian mask by the polynomial (eqn. 37)
v = bsxfun(@times, gauss_mask(x, y, r), polyval(coeff, v));
end