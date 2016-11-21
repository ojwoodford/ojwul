%DOG Difference of Gaussians blob detector
%
%   features = dog(I)
%
%IN:
%   I - HxWxC image
%
%OUT:
%   features - 4xN frames of N features: [X; Y; scale; orientation].

function features = dog(I)
if size(I, 3) == 3
    I = rgb2gray(I);
end
features = vl_sift(single(I));
end 
