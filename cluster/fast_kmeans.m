%FAST_KMEANS Compute k-means cluster centres
%
%	[CX sse I] = fast_kmeans(X, params[, CX])
%
% This function computes k-means cluster centres. It will work faster if
% data is projected onto principle components first.
% 
% IN:
% 	X - MxN matrix of N input vectors of dimension M.
%   params - [nclusters max_iters min_sse_delta robust_thresh].
%             nclusters: number of required clusters.
%             max_iters: maximum number of iterations.
%             min_sse_delta: convergence threshold for change in error score.
%             robust_thresh: reject outliers more than robust_thresh
%                            squared distance away from all cluster
%                            centers. If 0, no outlier rejection.
%   CX - Mx(nclusters) matrix of initial cluster centres of X. Default:
%        randomly initialized to nclusters different vectors from X.
% 
% OUT:
% 	CX - Mx(nclusters) matrix of updated cluster centres of X.
% 	sse - Sum of squared errors of distances to cluster centres.
% 	I - Nx1 uint32 vector of cluster indices each input vector is assigned
% 	    to.

function varargout = fast_kmeans(varargin)
sourceList = {'fast_kmeans.cpp', '-Xopenmp'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return