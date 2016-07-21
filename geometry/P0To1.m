%P0TO1 Convert projection matrices' principal points from 0 to 1 based
%
%   P = P0To1(P)

function P = P0To1(P)
for a = 1:size(P, 3)
    P(:,:,a) = [1 0 1; 0 1 1; 0 0 1] * P(:,:,a);
end