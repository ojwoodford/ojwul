%OJW_BOX_PYRAMID  Creates an image pyramid using box filtering
%
%   B = ojw_box_pyramid(A, max_depths)
%
% Generates the image pyramid of an image of any class.
%
% IN:
%   A - MxNxC input image, or Sx1 cell array of input images.
%   max_depths - scalar indicating maximum number of additional pyramid
%                levels to calculate.
%
% OUT:
%   B - 1xL (or SxL) cell array, where B{1} == A, and B{2},..,B{L} are the
%       subsampled pyramid levels.


function B = ojw_box_pyramid(A, max_depths)

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
else
    C = A;
end

% Calculate other levels
for a = 2:d
    % Filter and subsample
    if mod(size(C, 1), 2)
        C = padarray(C, [1 0], 'replicate', 'post');
    end
    C = C(1:2:end,:,:) + C(2:2:end,:,:);
    if mod(size(C, 2), 2)
        C = padarray(C, [0 1], 'replicate', 'post');
    end
    C = C(:,1:2:end,:) + C(:,2:2:end,:);
    C = C * 0.25;
    
    % Cast image
    if strcmp(in_class, 'logical')
        B{a} = C > 0.5;
    else
        B{a} = cast(C, in_class);
    end
end
return
