%COMPUTE_EDGELS  Extract edgels from a gradient image
%
%   E = compute_edgels(G, thresh)
%
% Returns a list of edgel positions and normal angles, each edgel being a
% pixel long, given a gradient image and threshold parameter.
%
% IN:
%   G - HxWx2 array of gradient images in x and y directions, along third
%       dimension.
%   thresh - Scalar value determining the threshold of the second
%            derivative of gradient at which to ignore edges.
%
% OUT:
%   E - 3xN list of edgels, each column giving [x y angle] of the edgel
%       centre.
%   I - Mx1 list of chain start indices (if min_length > 0)

function [E, I] = compute_edgels(G, thresh, min_length, join_thresh)
if nargin < 4
    join_thresh = 0;
    if nargin < 3
        min_length = 0;
        if nargin < 2
            thresh = 0;
        end
    end
end

% Compute the gradient maxima and subpixel locations
[M, E] = edge_grad_max(G, thresh);

% Compute the edgeness map
M_ = M(:,:,1) ~= 0;
M = (M(:,:,3) .* M_) ./ (M(:,:,2) + 20);

% Output only those edgelets which meet the criteria
G = reshape(G, [], 2);
E = [E(:,M_); normalize([-G(M_,2)'; G(M_,1)']); M(M_)'];

% Compute the edgelets or edge chains
I = [];
if min_length
    % Create an index matrix
    I = M;
    M_ = find(M_);
    n = numel(M_);
    I(M_) = 1:n;
    I(I==0) = n + 1;
    
    % Extract the neighbours to the right and below
    Y = [E zeros(5, 1)];
    off = size(M, 1);
    off = [-off-1, -off, -off+1, -1, 1, off-1, off, off+1]';
    I = I(bsxfun(@plus, M_', off));
    Y = reshape(Y(:,I), 5, 8, n);
    
    % Compute the edge score for the neighbours
    Z = normalize(bsxfun(@minus, reshape(E(1:2,:), 2, 1, []), Y(1:2,:,:)));
    Z = reshape(bsxfun(@times, abs(sum(bsxfun(@times, reshape(E(3:4,:), 2, 1, []), Z), 1) .* sum(Y(3:4,:,:) .* Z, 1)) .* Y(5,:,:), reshape(E(5,:), 1, 1, [])), 8, n);
    
    % Compute the chains
    [J, I] = compute_chains(Z, I, min_length, join_thresh);
    E = E(:,J);
end

% Compute the edge angles, in radians
E = [E(1:2,:); atan2(-E(3,:), E(4,:)); E(5,:)];
end

function [I, J] = compute_chains(Z, L, min_length, thresh)
I = zeros(size(Z, 2), 1);
ic = 0;
ic_last = 0;
J = I;
jc = 1;
K = uint8(sum(Z > thresh, 1)');
M = [K < 0; true];
i_ = uint32(0);
opts = {'==', uint8(1)};
while 1
    % Select the start point
    i = find_first(K, opts{:}, i_);
    if ~i
        if isempty(opts)
            break;
        end
        opts = {};
        continue;
    end
    i_ = i;
    J(jc) = ic + 1;
    % Build the chain
    l = 0;
    m = thresh * (i ~= uint32(0)) + 1;
    while m > thresh
        % Add to the chain
        ic = ic + 1;
        I(ic) = i;
        M(i) = true;
        K(i) = 0;
        l = l + 1;
        % Find the next edge
        k = Z(:,i);
        k(M(L(:,i))) = 0;
        [m, j] = max(k);
        i = L(j,i);
    end
    % Check it's long enough
    if l < min_length;
        ic = ic_last;
    else
        jc = jc + 1;
        ic_last = ic;
    end
end
I = I(1:ic);
J = J(1:jc-1);
end

