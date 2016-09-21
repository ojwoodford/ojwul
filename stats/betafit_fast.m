function phat = betafit_fast(x)
%BETAFIT_FAST Maximum-likelihood estimate of beta distribution parameters
%
%   BETAFIT_FAST(X) Returns the maximum likelihood estimates of the
%   parameters of the beta distribution given the data in the vector, X.
%   The data must be real, and in the range (0,1).
%
%   The function independently fits beta distributions to the columns of X,
%   or if X is a cell array, to the cells of X.

% (C) Copyright, Oliver Woodford 2015

% Compute starting parameters
if iscell(x)
    c = [reshape(cellfun(@(x) mean(log(x(:))), x), 1, []); ...
         reshape(cellfun(@(x) mean(log1p(-x(:))), x), 1, [])];
else
    c = [mean(log(x), 1); mean(log1p(-x), 1)];
end
phat = 1 - exp(c([2 1],:));
phat = bsxfun(@times, phat, 0.5 ./ (sum(phat) - 1));

% Do 5 Newton steps
for iter = 1:5
    s = sum(phat);
    J = bsxfun(@minus, psi(phat) - c, psi(s));
    a = psi(1, phat(1,:));
    b = psi(1, phat(2,:));
    s = psi(1, s);
    J = bsxfun(@times, J, 1 ./ (a .* s - a .* b + b .* s));
    phat = phat + [sum([b - s; s] .* J); sum([s; a - s] .* J)];
end
