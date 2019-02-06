%FIND_FIRST Fast, vectorized version of find(A, 1, 'first')
%
%	B = find_first(A)
%	B = find_first(A, start)
%	B = find_first(A, operator, value)
%	B = find_first(A, operator, value, start)
%
% Find the first element in each column vector of an array which meets a
% particular comparison criterion.
%
% Example:
%   find_first(rand(10, 10), '>', 0.8)
%
%IN:
%	A - MxN real numeric array.
%	operator - string of a comparison operator, e.g. '>'. Default: '~='.
%	value - 1xN or scalar array of values to use in the comparison.
%	        Default: 0.
%   start - 1xN or scalar uint32 array of offsets from the start of each
%           vector at which to start searching.
%
%OUT:
%	B - 1xN uint32 array of indices of first element in each column vector
%	    meeting the criterion, or 0 if none are found.

function varargout = find_first(varargin)
sourceList = {'find_first.cpp', '-Xopenmp'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return
