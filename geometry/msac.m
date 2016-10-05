%MSAC Robustly fit a model to data with the MSAC algorithm
%
%   [model, sqDists] = msac(X, fittingFunc, distFunc, minSamples, distThresh, maxTrials)
%
%IN:
%   X - MxN data matrix, with each column being one datum, to which we are
%       seeking to fit a model.
%   fittingFunc - Handle to a function that fits a model to minSamples data
%                 from X. It is assumed that the function is of the form:
%                    model = fittingFunc(X);
%                 The fitting function can return multiple models in a cell
%                 array. If this function cannot fit a model, it should an
%                 empty matrix.
%    distFunc - Handle to a function that computes the squared distances
%               for each datum in X, given a model. It is assumed that the
%               function is of the form:
%                  sqDist = distFunc(model, X);
%    minSamples - The minimum number of samples from X required by
%                 fittingFunc to fit a model.
%    distThresh - The squared distance threshold used to decide whether the
%                 point is an inlier or not.
%    maxTrials - Maximum number of iterations. Default: 200.
%
%OUT:
%   model - The best scoring model found.
%   sqDists - 1xN vector of squared distances for the best model found.

function [bestM, bestD, stats] = msac(X, fittingFunc, distFunc, minSamples, distThresh, maxTrials)
confidence = 0.99999; % Desired probability of choosing at least one sample free from outliers
confidence = log1p(-confidence);
if nargin < 6
    % Assume 30% inliers
    maxTrials = confidence / (log1p(-(0.3 ^ minSamples)) + 1e-300);
end
bestM = [];
bestD = [];
[~, nPoints] = size(X);
if nPoints < minSamples
    return
end

numTrials = 0;
numFailures = 0;
bestScore = Inf;
N = maxTrials;
while numTrials < N
    % Select a set of points to try
    ind = randperm(nPoints, minSamples);

    % Fit model(s) to this random selection of data points
    M = fittingFunc(X(:,ind));
    if isempty(M)
        numFailures = numFailures + 1;
        if (numTrials + numFailures) >= maxTrials
            break;
        end
        continue;
    end
    if ~iscell(M)
        M = {M};
    end

    for a = 1:numel(M)
        % Compute the distances for each datum
        D = distFunc(M{a}, X);

        % Compute a score
        score = sum(min(D, distThresh));

        % Compare the model to the current best
        if score < bestScore
            % Compute a new number of trials
            pNoOutliers = sum(D < distThresh);
            if pNoOutliers == 0
                continue;
            end
            pNoOutliers = (pNoOutliers / nPoints) ^ minSamples;
            N = min(maxTrials, confidence / (log1p(-pNoOutliers) + 1e-300));

            % Store the new best model
            bestScore = score;
            bestD = D;
            bestM = M{a};
        end
    end
    numTrials = numTrials + numel(M);
end
if nargout > 2
    stats = [numTrials numFailures maxTrials];
end
end
