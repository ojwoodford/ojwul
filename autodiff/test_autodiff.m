function test_autodiff()

% Load an image
A = imread('peppers.png');
sz = size(A);

% Create a target image
[y, x] = ndgrid(linspace(sz(1)/2-50, sz(1)/2+50, 100), linspace(sz(2)/2-50, sz(2)/2+50, 100));
ref = ojw_interp2(A, x, y);

% Create a small homography
Hgt = eye(3) + randn(3) * 0.0002;
X = Hgt * [x(:)'; y(:)'; ones(1, numel(x))];

% Create a cost function
H = eye(3);
for a = 1:100
    % The magic line!!!!
    H_ = autodiff(H);
    
    % Apply the homography
    Y = H_ * X;
    Y = proj(Y);
    
    % Sample the image
    tgt = ojw_interp2(A, Y(1,:), Y(2,:));
    tgt = reshape(tgt, size(ref));

    % Compute the residuals
    r = tgt - ref;
    
    % Read out the *automagically* computed gradient of the residuals
    J = reshape(grad(r), 9, []);
    r = double(r);
    
    % Compute the gauss-newton step
    step = reshape(J' \ r(:), 3, 3);
    
    % Visualization
    subplot(131)
    sc(ref);
    subplot(132);
    sc(r, 'diff');
    subplot(133);
    sc(double(tgt));
    drawnow;
    
    % Check for convergence
    if norm(step(:)) < 1e-10
        break;
    end
    
    % Apply the step
    H = H - step;
end
end
