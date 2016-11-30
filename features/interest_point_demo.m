%INTEREST_POINT_DEMO  Compute and visualize the interest points of an image
%
%   interest_point_demo
%   interest_point_demo(A, [scale, [thresh, [detector]]])
%
%IN:
%   A - HxWxC image. Default: Use peppers.png
%   scale - scale of the interest points to be detected (if applicable).
%           Default: 2.
%   thresh - threshold on point suppression. If negative, magnitude is
%            interpreted as the proportion of strongest points to keep.
%            Default: -0.5, i.e. keep top 50% of points.
%   detector - string of the detector name. Default: 'harris'.

function interest_point_demo(A, detector, scale, thresh)
% Set default arguments
if nargin < 4
    thresh = -0.5; % Negative numbers indicate ratio of edges to keep
    if nargin < 3
        scale = 2;
        if nargin < 2
            detector = 'harris';
            if nargin < 1
                A = imread('peppers.png');
            end
        end
    end
end

switch lower(detector)
    case 'harris'
        % Compute the detector score image
        tic;
        score = harris(A, scale);
        t1 = toc;
        
        % Compute the interest points
        tic;
        X = extract_features(score, 1.5, thresh, true);
        t2 = toc;
        
        % Visualization
        fprintf('Computing Harris detector score: %gs\nExtracting interest points: %gs\n', t1, t2);
        clf;
        % Render the score over the image
        sc(cat(3, score, A), 'prob');
        % Render the interest points
        hold on
        plot(X(1,:), X(2,:), 'y.', 'MarkerSize', 10);
        hold off
        
    case 'dog'
        % Compute the detector score image
        tic;
        features = dog(A);
        t1 = toc;
        
        % Visualization
        fprintf('Detecting DoG features: %gs\n', t1);
        clf;
        % Show the image
        imdisp(A);
        % Render the interest points
        hold on
        render_circles(features(1:3,:), 64, [1 0 0]);
        hold off
    otherwise
        error('Detector method %s not recognized', detector);
end
end
