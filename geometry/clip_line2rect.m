%CLIP_LINE2RECT Clips a line segment to within a rectangle
%
%   X = clip_line2rect(X, rect)
%
% Clip a line to a rectangle using the Cohen-Sutherland algorithm,
% described here: https://en.wikipedia.org/wiki/Cohen-Sutherland_algorithm
%
%IN:
% X - 2x2 matrix of [x1 x2; y1 y2] start and end points.
% rect - 2x2 matrix of [left right; top bottom] rectangle boundaries.
%
%OUT:
% X - 2x2 matrix of [x1 x1; y1 y2] clipped start and end points.

function X = clip_line2rect(X, rect)
c(1) = area_code(X(:,1), rect);
c(2) = area_code(X(:,2), rect);

while true
    if bitand(c(1), c(2))
        % Segment is out of the rectangle
        X(:) = NaN;
        return;
    end
    [c_, ind] = max(c);
    if  c_ == 0
        % Segment is inside the rectangle
        return;
    end
    % Clip to one of the crossing axes
    if c_ >= 8
        X(1,ind) = X(1,ind) + (X(1,3-ind) - X(1,ind)) * (rect(2,2) - X(2,ind)) / (X(2,3-ind) - X(2,ind));
        X(2,ind) = rect(2,2);
    elseif c_ >= 4
        X(1,ind) = X(1,ind) + (X(1,3-ind) - X(1,ind)) * (rect(2,1) - X(2,ind)) / (X(2,3-ind) - X(2,ind));
        X(2,ind) = rect(2,1);
    elseif c_ >= 2
        X(2,ind) = X(2,ind) + (X(2,3-ind) - X(2,ind)) * (rect(1,2) - X(1,ind)) / (X(1,3-ind) - X(1,ind));
        X(1,ind) = rect(1,2);
    else
        X(2,ind) = X(2,ind) + (X(2,3-ind) - X(2,ind)) * (rect(1,1) - X(1,ind)) / (X(1,3-ind) - X(1,ind));
        X(1,ind) = rect(1,1);
    end
    c(ind) = area_code(X(:,ind), rect);
end
end

function c = area_code(X, rect)
c = uint8(0);
if X(1) < rect(1,1)
    c = uint8(1);
elseif X(1) > rect(1,2)
    c = uint8(2);
end
if X(2) < rect(2,1)
    c = uint8(c) + 4;
elseif X(2) > rect(2,2)
    c = uint8(c) + 8;
end
end
    
