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
        d = d(1:end-1,:) / d(end,end);
        [~, d] = find_closest(Xt * d(:,1:end-1)' + d(:,end)');
        scores(a) = sum(d);
    end
    
    % Pick the best
    [~, a] = min(scores);
    T = T_ \ (R(:,:,a) * T);
    T = T(1:3,:) / T(4,4);
else
    % Start from the input
    T = eye(size(X, 1), size(X, 1)+1);
end
R = T(:,1:size(Y, 1))';
t = T(:,end);

% Debug rendering
if nargin < 4
    debug_vis = false;
elseif debug_vis
    figure();
    plot3(Yt(:,1), Yt(:,2), Yt(:,3), 'b.');
    hold on
    handle = plot3(Xt(:,1), Xt(:,2), Xt(:,3), 'r.');
end

% Iterate until convergence
best_score = Inf;
while true
    % Compute the closest points from X to Y
    X_ = Xt * R + t';
    closest = find_closest(X_);
    
    % Debug rendering
    if debug_vis
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
