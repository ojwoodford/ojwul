%HOMG Homogenize vectors
%
%   Y = homg(X)
%
% Homogenize an array of vectors along the first dimension, by
% concatenating a 1.
%
%IN:
%   X - MxN input array.
%
%OUT:
%   Y - (M+1)xN output array.

function X = homg(X)
X(end+1,:) = 1;