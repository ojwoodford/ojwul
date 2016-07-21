%P1TO0 Convert projection matrices' principal points from 1 to 0 based
%
%   P = P1To0(P)

function P = P1To0(P)
for a = 1:size(P, 3)
    P(:,:,a) = [1 0 -1; 0 1 -1; 0 0 1] * P(:,:,a);
end