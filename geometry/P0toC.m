%P0TOC Convert projection matrices' principal points from 0 to centre based
%
%   P = P0toC(P, im)

function P = P0toC(P, im)
if numel(im) == 2 || numel(im) == 3
    h = im(1);
    w = im(2);
else
    [h w c] = size(im);
end
K = [1 0 (1-w)/2; 0 1 (1-h)/2; 0 0 1];
for a = 1:size(P, 3)
    P(:,:,a) =  K * P(:,:,a);
end