%PCA Principal component analysis
%
%   [T, Y, V] = pca(X)
%
% IN:
%    X - MxN matrix of N vectors of dimension M.
%
% OUT:
%    T - Mx(M+1) Projection matrix for PCA transformation
%    Y - MxN matrix of transformed vectors, where Y = T * [X; ones(1, N)].
%    V - Mx1 list of eigen values associated with each dimension.

function [T, X, V] = pca(X)
% Check for NaNs
if ~any(isnan(X(:)))
    % Compute the mean and subtract
    t = mean(X, 2);
    X = bsxfun(@minus, X, t);

    % Do the economy svd, to avoid memory issues
    [R, V] = svd(X, 'econ');
    R = R';

    if nargout > 2
        % Square the singular values to give eigenvalues
        V = diag(V);
        V = (V .* V) ./ size(X, 2);
    end
else
    % Have NaNs!
    % Compute the mean and subtract
    t = mean(X, 2, 'omitnan');
    X = bsxfun(@minus, X, t);
    
    % Compute the covariance matrix
    sd = sqrt(mean(X .* X, 2, 'omitnan'));
    C = corrcoef(X', 'rows', 'pairwise');
    C = C .* (sd * sd');
    
    % Compute the eigen values and vectors
    [R, V] = eig(C);
    
    % Reorder so largest first
    V = flipud(diag(V));
    R = R(:,end:-1:1)';
end

% Compose the transform
T = [R -(R * t)];

if nargout > 1
    % Apply the transform
    X = R * X;
end
end