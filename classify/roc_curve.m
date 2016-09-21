%ROC_CURVE Compute the x and y parameters of an ROC curve
%
%   [X, Y, auc] = roc_curve(scores, ground_truth)
%
% Compute the parameters of a Receiver Operating Characteristic curve,
% which plots true positive rate against false positive rate, over a range
% of classification thresholds.
%
%IN:
%   scores - Mx1 vector of classification scores for M test variables.
%   ground_truth - Mx1 logical vector indicating the ground_truth class for
%                  each test variable (true is positive).
%
%OUT:
%   X - (M+1)x1 vector of false positive rates with increasing
%       classification threshold. 
%   Y - (M+1)x1 vector of true positive rates with increasing
%       classification threshold.
%   auc - scalar area under the curve.


function [X, Y, auc] = roc_curve(scores, ground_truth)

% Order the scores
[~, order] = sort(scores, 'descend');
ground_truth = ground_truth(order);

% Compute the true positive rate
Y = [0; cumsum(ground_truth(:))];
Y = Y / Y(end);

% Compute the false positive rate
X = [0; cumsum(~ground_truth(:))];
X = X / X(end);

if nargout < 3
    return
end

% Compute the area under the curve
auc = (diff(X)' * (Y(1:end-1) + Y(2:end))) * 0.5;
end