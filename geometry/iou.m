%IOU  Compute the intersection over union of two shapes
%
%   score = iou(A, B)
%
% Compute the intersection over union of two shapes, represented either by
% logical images, or by polygons.
%
%IN:
%   A - 2xM polygon coordinates, or HxW logical matrix.
%   B - 2xN polygon coordinates, or HxW logical matrix.
%
%OUT:
%   score - scalar intersection over union score, in the range [0, 1].

function score = iou(A, B)
if isnumeric(A) && size(A, 1) == 2 && isnumeric(B) && size(B, 1) == 2
    % Polygons
    A = polyshape(A');
    B = polyshape(B');
    score = area(intersect(A, B)) / area(union(A, B));
elseif islogical(A) && islogical(B) && isequal(size(A), size(B))
    % Bitmaps
    score = sum(A(:) & B(:)) / sum(A(:) | B(:));
else
    error('Inputs not recognized');
end
end
