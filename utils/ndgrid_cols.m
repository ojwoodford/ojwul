%NDGRID_COLS Like NDGRID, but creates column vectors from the outputs
%
%   [X, sz] = ndgrid_cols(...)
%
% This function applies passes its inputs directly to NDGRID, then converts
% the outputs to row vectors, which are stacked vertically, so each
% combination of inputs becomes a column vector in the output matrix.
%
%IN:
%   varargin - Same as inputs to NDGRID.
%
%OUT:
%   X - (nargin)xN matrix of permutation column vectors.
%   sz - 1xM array of sizes of outputs that ndgrid would produce with the
%        same inputs, such that N = prod(sz).

function [X, sz] = ndgrid_cols(varargin)
X = cell(nargin, 1);
[X{:}] = ndgrid(varargin{:});
sz = size(X{1});
X = cell2mat(cellfun(@(x) x(:)', X, 'UniformOutput', false));
end