%GAUSS_MASK  Compute 1D Nth derivative of a Gaussian
%
% Examples:
%    F = gauss_mask(sigma)
%    F = gauss_mask(sigma, deriv)
%    F = gauss_mask(sigma, deriv, X)
%
% This function computes the Nth derivative of a Gaussian, in 1D.
%
% IN:
%    sigma - Standard deviation of the Gaussian.
%    deriv - The derivative to compute. Default: 0.
%    X - Values of x at which to compute Gaussian. Default:
%        -floor(4*sigma):floor(4*sigma).
%
% OUT:
%    F - The vector of values of Nth derivative of the zero-mean Gaussian
%        with standard deviation sigma, evaluated at the values in X.

function F = gauss_mask(sigma, deriv, X)
if nargin < 3
    % Determine necessary filter support (for Gaussian).
    X = max(floor(4 * sigma), 1);
    X = -X:X;
end

% Compute standard Gaussian
s2 = -1 ./ sigma .^ 2;
F = exp((0.5 * s2) * X .* X) ./ (sigma * sqrt(2 * pi));

if nargin > 1 && deriv
    % Compute 1st derivative
    F_ = F;
    F = s2 * X .* F;
    % Compute nth derivative
    for a = 2:deriv
        F__ = F_;
        F_ = F;
        F = s2 * ((a - 1) .* F__ + X .* F_);
    end
end
return