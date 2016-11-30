%FAST_CORNERS  Call mexed FAST corner detector
%
%	XY = fast_corners(I, thresh, type)
%
% FAST corner detection using Ed Rosten's C implementation. Method
% published in:
%   "Machine learning for high-speed corner detection", 
%   E. Rosten & T. Drummond, ECCV 2006.
%
%IN:
%	I - HxW uint8 grayscale image.
%   thresh - scalar integer threshold for nonmax-suppression.
%   type - {9,10,11,12} type of corner detector to use.
%
%OUT:
%	XY - 2xN uint16 array of x,y coordinates of corners.

function varargout = fast_corners(varargin)
sd = 'private/fast/';
sourceList = {[sd 'fast_mex.c'], [sd 'fast_9.c'], [sd 'fast_10.c'], [sd 'fast_11.c'], ...
              [sd 'fast_12.c'], [sd 'nonmax.c'], [sd 'fast.c']}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
end