%PCTO0 Convert projection matrices' principal points from centre to 0 based
%
%   P = PcTo0(P, im)

function P = PcTo0(P, im)
if numel(im) == 2 || numel(im) == 3
    h = im(1);
    w = im(2);
else
    [h, w, c] = size(im);
end
K = [1 0 (w-1)/2; 0 1 (h-1)/2; 0 0 1];
P = tmult(K, P);
