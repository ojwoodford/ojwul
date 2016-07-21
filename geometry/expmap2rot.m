%EXPMAP2ROT Computes a rotation matrix from exponential map angle.
% 
% R = expmap2rot(X)

function R = expmap2rot(X)

% Compute the angle
sz = size(X);
X = reshape(X, 3, []);
angle = sqrt(sum(X .* X, 1));

% Initialize output
R = repmat([1 0 0 0 1 0 0 0 1]', 1, numel(angle));

% Only update rotations with non-zero angles
M = angle ~= 0;
angle = angle(M);

if ~isempty(angle)
    % Compute the axes, sin and cos
    axis = bsxfun(@times, X(:,M), 1./angle);
    
    % Compute the skew matrix
    cpm = zeros(9, numel(angle));
    cpm([2 6 7],:) = axis([3 1 2],:);
    cpm([3 4 8],:) = -axis([2 3 1],:);
    
    % Compute the squared skew matrix
    cpm2 = reshape(bsxfun(@times, reshape(axis, 3, 1, []), reshape(axis, 1, 3, [])), 9, []);
    cpm2([1 5 9],:) = bsxfun(@minus, cpm2([1 5 9],:), sum(cpm2([1 5 9],:)));
    
    % Apply Rodrigues' formula
    R(:,M) = R(:,M) + bsxfun(@times, cpm, sin(angle)) + bsxfun(@times, cpm2, 1 - cos(angle));
end

% Reshape for output
R = reshape(R, [3 3 sz(2:end)]);
