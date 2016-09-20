%VGG_ARGPARSE  Parse variable arguments into a structure
%  opts = vgg_argparse(inopts,varargin)
%    inopts: structure (cells array) listing valid members and default values
%    varargin: variable arguments of form '<name>',<value>,...
%    opts: opts modified by varargin
%
%  Example:
%    function f = foo(varargin)
%    opts = vgg_argparse(struct('maxiters',10,'verbose',0), varargin)
%    ...
%
%  An unknown option (ie, present in varargin but absent in inopts)
%  causes an error. Calling the function as 
%  [opts,rem_opts] = vgg_argparse(inopts,varargin) returns the unknown
%  option(s) in rem_opts for later use rather than causes an error.
%
%  May also use OPTS = VGG_ARGPARSE(OPTS, ASTRUCT) where ASTRUCT is a struct
%  of options.

% Author: Mark Everingham <me@robots.ox.ac.uk>
% modified by werner, Jan 03
% modified by ojw, Mar 14
% Date: 16 Jan 02

function [opts, rem_opts] = vgg_argparse(opts, varargin)

if iscell(opts)
  opts = struct(opts{:});
end

if isempty(varargin)
    inopts = struct([]);
elseif nargin == 2
    if isempty(varargin{1})
        inopts = struct([]);
    elseif isstruct(varargin{1})
        inopts = varargin{1};
    elseif iscell(varargin{1})
        inopts = struct(varargin{1}{:});
    else
        error('Single argument');
    end
else
    inopts = struct(varargin{:});
end

rem_opts = [];
for fn = fieldnames(inopts)'
    if isfield(opts, fn{1})
        opts.(fn{1}) = inopts.(fn{1});
    else
        rem_opts.(fn{1}) = inopts.(fn{1});
    end
end

if nargout < 2 && ~isempty(rem_opts)
    v = fieldnames(rem_opts);
    error(['Bad arguments: ' sprintf('''%s'' ', v{:})]);
end