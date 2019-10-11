% [R,s,t,Y1,p] = procrustes(X,Y) Procrustes alignment
%
% Finds the best similarity transformation Y = s.X.R + ones(N,1).t' (in the
% least squares sense).
%
% References:
% - Borg & Groenen: "Modern Multidimensional Scaling: Theory and Application",
%   Springer, 2005 (chapter 20).
% - Cox & Cox: "Multidimensional Scaling", 2nd ed. Chapman & Hall, 2000
%   (chapter 5).
%
% In:
%   X,Y: NxD data sets of row vectors.
% Out:
%   R: DxD orthogonal matrix.
%   s: scale.
%   t: Dx1 translation vector.
%   Y1: NxD matrix, transformed X.
%   p: Procrustes statistic (normalised error |Y-Y1|?/tr(Y.Y')).

% Copyright (c) 2009 by Miguel A. Carreira-Perpinan

function [R,s,t,Y1,p] = procrustes(X,Y)

[N, D] = size(Y);
sY = sum(Y,1); sX = sum(X,1); C = Y'*X - sY'*(sX/N); [U,S,V] = svd(C);
R = V*U';
s = sum(sum(C.*R')) / (sum(sum(X.^2))-(sX*sX')/N);
t = (sY - s*sX*R)'/N;

if nargout > 3
  Y1 = bsxfun(@plus,s*X*R,t');
  p = sum(sum((Y1-Y).^2))/sum(sum(Y.^2));
end
end