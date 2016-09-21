%FIT_GAUSSIAN Fit a multi-variate gaussian to data.
%
%   [mu, whiten] = fit_gaussian(X)
%
% Fit a multi-variate gaussian to data. This function can handle
% under-constrained data.
%
% IN:
%    X - MxN matrix of N vectors of dimension M.
%
% OUT:
%    mu - Mx1 distribution mean.
%    whiten - Mx(min(M,N)) whitening matrix.

function [mu, whiten] = fit_gaussian(X)
% Compute the mean and subtract
mu = mean(X, 2);
X = bsxfun(@minus, X, mu);

% Do the economy svd, to avoid memory issues
[whiten, V] = svd(X, 'econ');

% Compute the whitening matrix
whiten = diag(sqrt(size(X, 2) - 1) ./ (diag(V) + 1e-300)) * whiten';
end