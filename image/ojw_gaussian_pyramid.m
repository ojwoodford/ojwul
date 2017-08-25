%OJW_GAUSSIAN_PYRAMID  Creates a gaussian image pyramid
%
%   B = ojw_gaussian_pyramid(A, max_depths, [F])
%
% Generates the gaussian pyramid of an image of any class.
%
% IN:
%   A - MxNxC input image, or Sx1 cell array of input images.
%   max_depths - scalar indicating maximum number of additional pyramid
%                levels to calculate.
%   F - Tx1 filter to apply to image in x and y directions before
%       downsampling. Default: fspecial('gaussian', [5 1], 1/sqrt(2)).
%
% OUT:
%   B - 1xL (or SxL) cell array, where B{1} == A, and B{2},..,B{L} are the
%       subsampled pyramid levels.


function B = ojw_gaussian_pyramid(A, max_depths, F)

if nargin < 3
    % Use a separable smoothing filter
    F = [0.25 0.5 0.25];
end

% Multi-cell input
if iscell(A)
    B = reshape(A, numel(A), 1);
    B{1,max_depths+1} = [];
    for a = 1:numel(A)
        C = ojw_gaussian_pyramid(A{a}, max_depths, F);
        B(a,1:numel(C)) = C;
    end
    return
end

% Calculate number of depths required
d = ceil(log2(max(size(A, 1), size(A, 2))));
if nargin > 1
    d = min(d, max_depths);
end
d = d + 1;

% Initialize the first level
B = cell(1, d);
B{1} = A;

% Convert to singles for higher accuracy (especially with logical arrays)
% but leave doubles as they are
in_class = class(A);
if ~isa(A, 'double')
    C = single(A);
    F = single(F);
else
    C = A;
    F = double(F);
end

% Calculate other levels
for a = 2:d
    % Filter and subsample
    C = filter_subsample(C, F);
    
    % Cast image
    if strcmp(in_class, 'logical')
        B{a} = C > 0.5;
    else
        B{a} = cast(C, in_class);
    end
end
return
