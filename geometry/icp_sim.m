%ICP_SIM  Compute the similarity transform that best aligns two point sets
%
%    T = icp_sim(X, Y, [initialize])
%
% Given two sets of points, X and Y, this function iteratively solves the
% optimization problem:
%
%  T = argmin_T sum_i min_j || T * homg(X(:,i)) - Y(:,j) || ^ 2
%
% subject to T being a similarity transform.
%
%IN:
%   X - DxM set of D-dimensional points to transform
%   Y - DxN set of D-dimensional points to align to (which should be well
%       distributed)
%   initialize - Scalar, 0 (no initialization), 1 (try 4 starting
%                positions), 2 (optimize 4 starting positions for 4
%                iterations each).
%
%OUT:
%   T - Dx(D+1) similarity transform matrix.

function T = icp_sim(X, Y, initialize, debug_vis)
% Initialize the search function
Xt = X';
Yt = Y';
tri = delaunayn(Yt);
find_closest = @(X) dsearchn(Yt, tri, X);

% Initialize the output
if nargin > 2 && initialize
    % Prealign using principle components
    [~, T] = whiten_srt(X);
    [~, T_] = whiten_srt(Y);
    
    % Try 4 possible axis flips
    R = cat(3, diag([1 1 1 1]), diag([1 -1 -1 1]), diag([-1 1 -1 1]), diag([-1 -1 1 1]));
    for a = 4:-1:1
        d = T_ \ (R(:,:,a) * T);
        T__{a} = d(1:end-1,:) / d(end,end);
        if initialize == 1
            [~, d] = find_closest(Xt * T__{a}(:,1:end-1)' + T__{a}(:,end)');
            scores(a) = sum(d);
        else
            [T__{a}, scores(a)] = run_iters(Xt, Yt, T__{a}, find_closest, [], 4);
        end
    end
    
    % Pick the best
    [~, a] = min(scores);
    T = T__{a};
else
    % Start from the input
    T = eye(size(X, 1), size(X, 1)+1);
end

% Debug rendering
handle = [];
if nargin > 3 && debug_vis
    figure();
    plot3(Yt(:,1), Yt(:,2), Yt(:,3), 'b.');
    hold on
    handle = plot3(Xt(:,1), Xt(:,2), Xt(:,3), 'r.');
end

% Iterate until convergence
T = run_iters(Xt, Yt, T, find_closest, handle, 100);
end

function [T, best_score] = run_iters(Xt, Yt, T, find_closest, handle, iters)
R = T(:,1:end-1)';
t = T(:,end);
best_score = Inf;
while iters > 0
    iters = iters - 1;
    
    % Compute the closest points from X to Y
    X_ = Xt * R + t';
    closest = find_closest(X_);
    
    % Debug rendering
    if ~isempty(handle)
        set(handle, 'XData', X_(:,1), 'YData', X_(:,2), 'ZData', X_(:,3));
        drawnow();
    end
    
    % Solve the procrustes problem
    [R, s, t, ~, score] = procrustes(Xt, Yt(closest,:));
    
    % Check the score decreased
    if score >= best_score
        break;
    end
    best_score = score;
    
    % Compute the transform
    R = s * R;
    T = [R', t];
end
end
