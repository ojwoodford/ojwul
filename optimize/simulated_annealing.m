%SIMULATED_ANNEALING  Perform simulated annealing on conditional energies
%
%   [X en] = simulated_annealing(X0, energy, T, varargin)
%
% Simulated annealing for discrete energy minimization labelling problems
% for which conditional energy distributions can be computed. Given an
% initial labelling, a normalized cooling schedule, and a function which
% returns conditional energy distributions, perform simulated annealing,
% cycling over all variables in a fixed, random order, once for each input
% temperature.
%
% The cooling schedule is normalized so that 100 means random (uniform)
% sampling, while 0 (freezing) means deterministic selection of the energy
% minimizing label (i.e. ICM).
%
%IN:
%   X0 - Initial labelling array, containing for each variable in the state
%        space a label index from 1 to n, where n is the number of possible
%        labels for that variable (a value which can vary from variable to
%        variable).
%   energy - A handle to a function to compute the current total or
%            conditional energy, with the following arguments:
%                  E = energy(X, j, varargin{:});
%            X is given labelling, j is the index of the variable in X for
%            which the conditional energy distribution, E, is to be
%            computed (E has length n, where n is the number of possible
%            labels for variable j). If j = 0 then the function returns the
%            energy of the state X (so E is scalar). varargin are the
%            trailing input variables passed into the function
%            simulated_annealing.
%   T - Nx1 normalized temperature cooling schedule, with a temperature
%       value for each cycle through all of the input variables (i.e. N
%       cycles total).
%   varargin - All trailing input variables are passed to the function
%              energy. Anonymous functions can be used instead, but this
%              method can be faster.
%
%OUT:
%   X - Lowest energy labelling found.
%   en - Energy of X.

% $Id: simulated_annealing.m,v 1.9 2008/10/03 15:59:04 ojw Exp $

function [X en] = simulated_annealing(X0, energy, T, varargin)

% Calculate initial energy
en = zeros(numel(T)+1, 1);
en(1) = energy(X0, 0, varargin{:});
enMin = en(1);
X = X0;

% Cache values needed to evaluate start temperature
nX = numel(X0);
nT = numel(T);
if any(T(:))
    Tcurr = -1;
    Tot = 0;
    sumN = 0;
    for a = 1:nX
        % Calculate the conditional energy
        enCon = energy(X0, a, varargin{:});
        Tot = Tot + sum(exp(min(enCon) - enCon));
        sumN = sumN + numel(enCon);
    end
    % Scale temperatures
    T = (min(max(T, 0), 100) * (0.01 * (sumN - nX)) + nX) / sumN;
    T0 = nX / sumN;
else
    % Just doing ICM
    T0 = 0;
end

I = randperm(nX); % Use a random order, but the same one every iteration
for t = 1:nT
    % Set the temperature automatically
    if T(t) == T0
        % Freeze
        Tcurr = -Inf;
    else
        %fprintf('Ratio wanted: %g. Ratio got: %g.\n', T(t), Tot/sumN);
        % Estimate correct temperature for next iteration
        Tcurr = Tcurr * (log(T(t)) / log(Tot/sumN-eps));
        if Tcurr == -Inf
            % Stuck freezing when we shouldn't be
            Tcurr = -1; % Arbitrary!
        elseif Tcurr >= 0
            % Stuck trying to overcook
            Tcurr = -eps;
        end
    end
    en(t+1) = en(t);
    Tot = 0;
    for a = I
        % Calculate the conditional energy
        enCon = energy(X0, a, varargin{:});
        % Sample from the conditional distribution
        if Tcurr == -Inf
            % Freeze
            [s P] = min(enCon);
        else
            % Draw from the distribution, raised by some temperature
            P = enCon * Tcurr;
            P = cumsum(exp(P - max(P)));
            Tot = Tot + P(end);
            P = 1 + sum(P(1:end-1) < (rand(1) * P(end)));
        end
        en(t+1) = en(t+1) + (enCon(P) - enCon(X0(a)));
        X0(a) = P;
    end
    % Keep track of the lowest energy
    if en(t+1) < enMin
        X = X0;
        enMin = en(t+1);
    end
    ojw_progressbar('Simulated annealing', t / nT);
end
return