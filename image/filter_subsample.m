%FILTER_SUBSAMPLE  Filters and subsamples an image by a factor of 2
%
%   B = filter_subsample(A, F)
%
% Given an image and a one dimensional filter, this function convolves the
% image with the filter along its first two dimensions, and then subsamples
% the resulting image by a factor of two in both those dimensions.
%
% Support for new types, and for two filters for non-symmetrical, separable
% filtering can easily be added to this function.
%
% IN:
%   A - MxNxC input image.
%   F - Lx1 filter to be applied separably to A. F must be double iff A is
%       double, otherwise it must be single.
%
% OUT:
%   B - (M+1)/2x(N+1)/2xC output image, of the same class as A.

% $Id: filter_subsample.m,v 1.2 2007/08/02 11:34:56 ojw Exp $

function varargout = filter_subsample(varargin)
sourceList = {'filter_subsample.cpp', '-Xopenmp'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return