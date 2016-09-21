%DP_PAIR_CHAIN  Dynamic programming on a chain of pairwise links
%
% Examples:
%	L = dp_pair_chain(U, E)
%	[L en] = dp_pair_chain(U, E)
%
% Minimize the cost of a set of pairwise chains.
%
%IN:
%	U - PxQxR array of unary costs (double, single or uint32)
%	E - PxP or PxPx(Q-1)xR array of pairwise costs (same type as U).
%
%OUT:
%	L - QxR int32 array of optimal labels (numbered between 1 and P).
%   en - 1xR vector of minimum costs per chain.

% $Id: dp_pair_chain.m,v 1.1 2007/12/07 11:27:55 ojw Exp $

function varargout = dp_pair_chain(varargin)
sourceList = {'dp_pair_chain.cpp', '-Xopenmp'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return