%EXPM_SRT_3D Compute a transformation matrix, given the Lie algebra vector
%
%   M = expm_srt_3d(X)
%
% Computes the transformation matrix defined by a Lie vector consisting of
% a rotation, translation and uniform scaling.
%
% The computation is done in closed form, using the formulae given in the
% paper:
% "Distances and Means of Direct Similarities"
% M-T Pham et al.
%
%IN:
%   X - DxN array, where each column specifies a different
%       transformation. X(1:3,:) are the rotation components, X(4:6,:) are
%       the translation components and X(7,:) are the scale components. D
%       can be 3 (rotation only), 6 (rotation and translation), or 7
%       (rotation, translation and scale).
%
%OUT:
%   M - 3x(3+D~=3)xN array of transformation matrices.

function varargout = expm_srt_3d(varargin)
sourceList = {'expm_srt_3d.cpp', '-Xopenmp'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
