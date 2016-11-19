%IMFILTSEP Fast separable filtering on a multi-channel image
%
%   B = imfiltsep(A, fy, fx)
%
% Fast separable filtering of multi-channel images, using symmetric
% padding, providing an output of the same size as the input.
%
%IN:
%   A - HxWxC input image.
%   fy - filter to apply to the columns.
%   fx - filter to apply to the rows.
%
%OUT:
%   B - HxWxC filtered output image.

function I = imfiltsep(I, fy, fx)
% Compute the padding indices
[H, W, C] = size(I);
sympadding = @(N, n) [floor(n/2)+1:-1:2 1:N N-1:-1:N-ceil(n/2-0.5)];
Y = sympadding(H, numel(fy));
X = sympadding(W, numel(fx));

% For each channel separately
for c = 1:C
    % Compute the padded array
    J = I(Y,X,c);
    % Convolve
    J = conv2(fy(:), fx(:), J, 'valid');
    % Insert into array
    I(:,:,c) = J;
end
end
