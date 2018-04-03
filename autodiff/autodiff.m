% A class for doing autodifferentiation
classdef autodiff
    properties (SetAccess = private, Hidden = true)
        value;   % Function values    Mx...xN
        deriv;   % Jacobian values  VxMx...xN
        varind;  % Variable indices   1xV
    end
    
    methods
        function obj = autodiff(a, v, b)
            obj.value = a;
            if nargin > 1
                obj.varind = v(:)';
                %assert(isequal(obj.varind, unique(obj.varind)));
                if nargin > 2
                    obj.deriv = b;
                    %assert(isequal(size(b), [numel(obj.varind) size(a)]));
                else
                    % We assume that there is 1 varind per element in a
                    obj.deriv = reshape(eye(numel(v)), [numel(a) size(a)]);
                end
            else
                obj.varind = 1:numel(a);
                obj.deriv = reshape(eye(numel(a)), [numel(a) size(a)]);
            end
        end
        
        % Get the values
        function c = double(obj)
            c = obj.value;
        end
        
        function c = grad(obj, vars)
            if nargin < 2
                vars = obj.varind;
            elseif isscalar(vars) && vars < 0
                % Return sparse output
                n = numel(obj.value);
                c = sparse(repmat(obj.varind(:), [n 1]), ...
                           reshape(repmat(1:n, [numel(obj.varind(:)) 1]), [], 1), ...
                           obj.deriv(:), ...
                           -vars, n);
                return;
            end
            n = numel(vars);
            if n == numel(obj.varind) && n == vars(end) && n == obj.varind(end)
                c = obj.deriv;
            else
                c = zeros([n size(obj.value)]);
                [Lia, Locb] = ismember(vars, obj.varind);
                c(Lia,:) = obj.deriv(Locb(Lia),:);
            end
        end
        
        function c = var_indices(obj)
            c = obj.varind;
        end
        
        function disp(obj, varargin)
            disp(obj.value, varargin{:});
        end
        
        % Elementwise operators
        function c = plus(a, b)
            c = bsxfun(@plus, double(a), double(b));
            [d, v] = combine_grads(grad(a), grad(b), var_indices(a), var_indices(b), size(c), 'add');
            c = autodiff(c, v, d);
        end
        
        function c = uminus(a)
            c = autodiff(-a.value, a.varind, -a.deriv);
        end
        
        function c = minus(a, b)
            c = a + (-b);
        end
        
        function c = times(a, b)
            da = double(a);
            db = double(b);
            ga = [];
            gb = [];
            if isautodiff(a)
                ga = bsxfun(@times, a.deriv, shiftdim(db, -1));
            end
            if isautodiff(b)
                gb = bsxfun(@times, b.deriv, shiftdim(da, -1));
            end
            c = bsxfun(@times, da, db);
            [d, v] = combine_grads(ga, gb, var_indices(a), var_indices(b), size(c), 'add');
            c = autodiff(c, v, d);
        end
        
        function c = rdivide(a, b)
            da = double(a);
            db = double(b);
            ga = [];
            gb = [];
            if isautodiff(a)
                ga = bsxfun(@times, a.deriv, shiftdim(db, -1));
            end
            if isautodiff(b)
                gb = bsxfun(@times, b.deriv, -shiftdim(da, -1));
            end
            c = bsxfun(@rdivide, da, db);
            [d, v] = combine_grads(ga, gb, var_indices(a), var_indices(b), size(c), 'add');
            c = autodiff(c, v, bsxfun(@times, d, shiftdim(1 ./ (db .* db), -1)));
        end
        
        function c = ldivide(a, b)
            c = rdivide(b, a);
        end
        
        function c = bsxfun(func, a, b)
            c = func(a, b);
        end
        
        function c = conj(a)
            c = autodiff(conj(a.value), a.varind, conj(a.deriv));
        end
        
        function c = abs(a)
            M = a.value < 0;
            c = reshape(a.deriv, [], numel(a.value));
            c(:,M) = -c(:,M);            
            c = autodiff(abs(a.value), a.varind, reshape(c, size(a.deriv)));
        end
         
        function c = exp(a)
            c = exp(a.value);
            c = autodiff(c, a.varind, bsxfun(@times, shiftdim(c, -1), a.deriv));
        end
        
        function c = log(a)
            c = autodiff(log(a.value), a.varind, bsxfun(@times, a.deriv, 1 ./ shiftdim(a.value, -1)));
        end
        
        function c = sqrt(a)
            c = sqrt(a.value);
            c = autodiff(c, a.varind, bsxfun(@times, a.deriv, shiftdim(0.5 ./ c, -1)));
        end
        
        function c = sin(a)
            c = autodiff(sin(a.value), a.varind, bsxfun(@times, a.deriv, shiftdim(cos(a.value), -1)));
        end
        
        function c = cos(a)
            c = autodiff(cos(a.value), a.varind, bsxfun(@times, a.deriv, shiftdim(-sin(a.value), -1)));
        end
        
        function c = tan(a)
            c = sec(a.value);
            c = autodiff(tan(a.value), a.varind, bsxfun(@times, a.deriv, shiftdim(c .* c, -1)));
        end
        
        function c = asin(a)
            c = autodiff(asin(a.value), a.varind, bsxfun(@times, a.deriv, shiftdim(1 ./ sqrt(1 - a.value .* a.value), -1)));
        end
        
        function c = acos(a)
            c = autodiff(acos(a.value), a.varind, bsxfun(@times, a.deriv, shiftdim(-1 ./ sqrt(1 - a.value .* a.value), -1)));
        end
        
        function c = atan(a)
            c = autodiff(atan(a.value), a.varind, bsxfun(@times, a.deriv, shiftdim(1 ./ (1 + a.value .* a.value), -1)));
        end
        
        % Matrix operators
        function c = mtimes(a, b)
            if isscalar(a) || isscalar(b)
                c = times(a, b);
                return;
            end
            da = double(a);
            db = double(b);
            c = da * db;
            v = unique([var_indices(a) var_indices(b)]);
            sz = [numel(v) size(c)];
            if isautodiff(a)
                g = reshape(sum(bsxfun(@times, grad(a, v), shiftdim(db, -2)), 3), sz);
            else
                g = 0;
            end
            if isautodiff(b)
                g = g + reshape(sum(bsxfun(@times, shiftdim(da, -1), permute(grad(b, v), [1 4 2 3])), 3), sz);
            end
            c = autodiff(c, v, g);
        end
        
        function c = mrdivide(a, b)
            if isscalar(a) || isscalar(b)
                c = rdivide(a, b);
                return;
            end
            error('Matrix divides not yet supported in autodiff. Use inv if possible.');
        end
        
        function c = mldivide(a, b)
            if isscalar(a) || isscalar(b)
                c = ldivide(a, b);
                return;
            end
            error('Matrix divides not yet supported in autodiff. Use inv if possible.');
        end
        
        function c = inv(a)
            c = inv(a.value);
            sz = [size(a.deriv) 1];
            d = bsxfun(@times, -shiftdim(c, -1), reshape(a.deriv, sz([1 end 2:end-1])));
            d = reshape(sum(d, 3), sz);
            d = bsxfun(@times, d, shiftdim(c, -2));
            d = reshape(sum(d, 3), sz);
            c = autodiff(c, a.varind, d);
        end
        
        function c = expm(a)
            c = expm(a.value);
            sz = [size(a.deriv) 1];
            d = bsxfun(@times, shiftdim(c, -1), reshape(a.deriv, sz([1 end 2:end-1])));
            d = reshape(sum(d, 3), sz);
            c = autodiff(c, a.varind, d);
        end
        
        function c = logm(a)
            sz = [size(a.deriv) 1];
            d = bsxfun(@times, shiftdim(inv(a.value), -1), reshape(a.deriv, sz([1 end 2:end-1])));
            d = reshape(sum(d, 3), sz);
            c = autodiff(logm(a.value), a.varind, d);
        end
        
        % Reduction methods: Sum, prod, min, max
        function c = sum(a, dim, flag)
            assert(nargin < 3 || strcmp(flag, 'default'), 'Only default option supported');
            if nargin < 2
                dim = first_nonsingleton_dim(a);
            end
            c = autodiff(sum(a.value, dim), a.varind, sum(a.deriv, dim+1));
        end
        
        function c = prod(a, dim)
            if nargin < 2
                dim = first_nonsingleton_dim(a);
            end
            c = a.value;
            c(c==0) = 1e-154;
            c = autodiff(prod(a.value, dim), a.varind, sum(bsxfun(@times, a.deriv, shiftdim(bsxfun(@times, prod(c, dim), 1 ./ c), -1)), dim+1));
        end
        
        function c = select(a, b, M)
            da = double(a);
            db = double(b);
            v = unique([var_indices(a) var_indices(b)]);
            ga = grad(a, v);
            gb = grad(b, v);
            if isscalar(a)
                M = ~M;
                db(M) = da;
                gb(:,M) = repmat(ga, [1 sum(M(:))]);
                c = autodiff(db, v, gb);
            elseif isscalar(b)
                da(M) = db;
                ga(:,M) = repmat(gb, [1 sum(M(:))]);
                c = autodiff(da, v, ga);
            else
                da(M) = db(M);
                ga(:,M) = gb(:,M);
                c = autodiff(da, v, ga);
            end
        end
        
        function [c, d] = min(a, b, dim)
            if nargin > 1 && ~isempty(b)
                d = double(a) > double(b);
                c = select(a, b, d);
            else
                if nargin < 3
                    dim = first_nonsingleton_dim(a);
                end
                [c, d] = min(a.value, [], dim);
                c = autodiff(c, var_indices(a), dimsel(a.deriv, d, dim));
            end
        end
        
        function [c, d] = max(a, b, dim)
            if nargin > 1 && ~isempty(b)
                d = double(a) < double(b);
                c = select(a, b, d);
            else
                if nargin < 3
                    dim = first_nonsingleton_dim(a);
                end
                [c, d] = max(a.value, [], dim);
                c = autodiff(c, var_indices(a), dimsel(a.deriv, d, dim));
            end
        end
        
        % Logical operators
        function c = lt(a, b)
            c = lt(double(a), double(b));
        end
        
        function c = le(a, b)
            c = le(double(a), double(b));
        end
        
        function c = gt(a, b)
            c = gt(double(a), double(b));
        end
        
        function c = ge(a, b)
            c = ge(double(a), double(b));
        end
        
        function c = eq(a, b)
            c = eq(double(a), double(b));
        end
        
        function c = ne(a, b)
            c = ne(double(a), double(b));
        end
        
        % Access functions
        function c = end(a, k, n)
            c = size(a.value, k);
        end
         
        function c = subsref(a, s)
            assert(strcmp(s.type, '()'));
            c = a.value(s.subs{:}); 
            c = autodiff(c, a.varind, reshape(a.deriv(:,s.subs{:}), [size(a.deriv, 1) size(c)]));
        end
        
        function c = subsasgn(a, s, b)
            assert(strcmp(s.type, '()'));
            c = double(a);
            c(s.subs{:}) = double(b);
            v = unique([var_indices(a) var_indices(b)]);
            d = grad(a, v);
            s_ = [{':'} s.subs];
            d_ = grad(b, v);
            if isscalar(b)
                d_ = repmat(d_, [1 size(c(s.subs{:}))]);
            end
            d(s_{:}) = d_;
            c = autodiff(c, v, d);
        end
        
        % Array concatenation
        function c = cat(dim, varargin)
            c = cellfun(@double, varargin, 'Uniform', false);
            v = cellfun(@var_indices, varargin, 'Uniform', false);
            v = unique(cat(2, v{:}));
            d = cellfun(@(x) grad(x, v), varargin, 'Uniform', false);
            c = autodiff(cat(dim, c{:}), v, cat(dim+1, d{:}));
        end
        
        function c = horzcat(varargin)
            c = cat(2, varargin{:});
        end
        
        function c = vertcat(varargin)
            c = cat(1, varargin{:});
        end
        
        % Transpose, permute, reshape, shiftdim etc.
        function c = ctranspose(a)
            c = permute(conj(a), [2 1]);
        end  
        
        function c = transpose(a)
            c = permute(a, [2 1]);
        end
        
        function c = permute(a, order)
            c = autodiff(permute(a.value, order), a.varind, permute(a.deriv, [1 order+1]));
        end
        
        function c = ipermute(a, order)
            order(order) = 1:numel(order);
            c = permute(a, order);
        end
        
        function c = reshape(a, varargin)
            if nargin == 2 && isnumeric(varargin{1})
                c = autodiff(reshape(a.value, varargin{1}), a.varind, reshape(a.deriv, [size(a.deriv, 1) varargin{1}]));
            else
                c = autodiff(reshape(a.value, varargin{:}), a.varind, reshape(a.deriv, size(a.deriv, 1), varargin{:}));
            end
        end
        
        function c = shiftdim(a, n)
            if (n > 0)
                m = ndims(a.value);
                n = rem(n, m);
                c = permute(a, [n+1:m, 1:n]);
            else
                c = reshape(a, [ones(1, -n), size(a.value)]);
            end
        end
        
        % Size and is*
        function varargout = size(a, varargin)
            [varargout{1:nargout}] = size(a.value, varargin{:});
        end
        
        function c = length(a)
            if isempty(a)
                c = 0;
            else
                c = max(size(a));
            end
        end
        
        function c = ndims(a)
            c = ndims(a.value);
        end
        
        function c = numel(a)
            c = numel(a.value);
        end
        
        function c = isautodiff(a)
            c = true;
        end
        
        function c = isempty(a)
            c = isempty(a.value);
        end
        
        function c = isscalar(a)
            c = isscalar(a.value);
        end
        
        function c = isreal(a)
            c = isreal(a.value) & isreal(a.deriv);
        end
        
        % Other functions        
        function c = ojw_interp2(I, x, y, interp_mode, oobv)
            if nargin < 5
                interp_mode = 'l';
            end
            if nargin < 5
                oobv = NaN;
            end
            assert(lower(interp_mode(1)) ~= 'c', 'Cubic interpolation not supported');
            if isautodiff(x) && isautodiff(y) && ~isautodiff(I)
                [c, d] = ojw_interp2(I, x.value, y.value, interp_mode, oobv);
                c = autodiff(c, x.varind, bsxfun(@times, x.deriv, d(1,:,:,:)) + bsxfun(@times, y.deriv, d(2,:,:,:)));
            elseif isautodiff(I) && ~isautodiff(x) && ~isautodiff(y)
                [h, w, n] = size(I.value);
                sz = size(I.value);
                [sz(1), sz(2)] = size(x);
                c = [shiftdim(I.value, -1); I.grad];
                c = reshape(reshape(c, size(c, 1), numel(I.value))', h, w, []);
                c = ojw_interp2(c, x, y, interp_mode, oobv);
                c = reshape(c, size(c, 1), size(c, 2), n, numel(I.varind)+1);
                c = autodiff(reshape(c(:,:,:,1), sz), I.varind, reshape(permute(c(:,:,:,2:end), [4 1 2 3]), [numel(I.varind) sz]));
            else
                error('Unexpected variables');
            end
        end
        
        function c = conv2(varargin)
            % Get the shape
            if ischar(varargin{end})
                shape = varargin{end};
                varargin = varargin(1:end-1);
            else
                shape = 'full';
            end
            v = combine_varind(varargin{:});
            n = numel(v);
            c = double(varargin{end});
            if isautodiff(varargin{end})
                d = permute(grad(varargin{end}, v), [2 3 1]);
            else
                d = [];
            end
            if numel(varargin) == 3
                varargin{1} = reshape(varargin{1}, [], 1);
                varargin{2} = reshape(varargin{2}, 1, []);
                conv2_ = @(A, B) conv2(B, A, shape);
            else
                conv2_ = @(A, B) conv2(A, B, shape);
            end
            for a = 1:numel(varargin)-1
                ca = double(varargin{a});
                if ~isempty(d)
                    for b = n:-1:1
                        d_(:,:,b) = conv2_(ca, d(:,:,b));
                    end
                    d = d_;
                    clear d_
                end
                if isautodiff(varargin{a})
                    da = permute(grad(varargin{a}, v), [2 3 1]);
                    for b = n:-1:1
                        d_(:,:,b) = conv2_(da(:,:,b), c);
                    end
                    if isempty(d)
                        d = d_;
                    else
                        d = d + d_;
                    end
                    clear d_
                end
                c = conv2_(ca, c);
            end
            d = permute(d, [3 1 2]);
            c = autodiff(c, v, d);
        end
        
        % Debug
        function check_size(a)
            sz1 = size(a.value);
            sz2 = size(a.deriv);
            if numel(sz2) == 2
                sz2(3) = 1;
            end
            assert(isequal(sz1, sz2(2:end)), 'Unexpected array sizes. Value: [%s]. Derivative: [%s].\n', sprintf('%d ', sz1), sprintf('%d ', sz2));
        end
    end
end

% Helpers
function dim = first_nonsingleton_dim(a)
dim = find(size(a) > 1, 1, 'first');
if isempty(dim)
    dim = 1;
end
end

function g = expand_array(ga, va, v, fill)
[l, l] = ismember(va, v);
sz = size(ga);
sz(1) = numel(v);
g = repmat(fill, sz);
g(l,:) = ga(:,:);
end

function [c, v] = combine_grads(ga, gb, va, vb, sz, M)
if ischar(M)
    switch M
        case 'add'
            func = @plus;
            fill = 0;
        case 'mult'
            func = @times;
            fill = 1;
    end
else
    func = @plus;
    fill = 0;
end
if isempty(va)
    sz2 = size(gb);
    sz = [sz2(1) sz];
    sz2(end+1:numel(sz)) = 1;
    c = repmat(gb, sz ./ sz2);
    v = vb;
    return;
end
if isempty(vb)
    sz2 = size(ga);
    sz = [sz2(1) sz];
    sz2(end+1:numel(sz)) = 1;
    c = repmat(ga, sz ./ sz2);
    v = va;
    return;
end
n = numel(va);
if n == numel(vb) && ((n == va(end) && n == vb(end)) || isequal(va, vb))
    c = bsxfun(func, ga, gb);
    v = va;
    return;
end
v = unique([va vb]);
c = bsxfun(func, expand_array(ga, va, v, fill), expand_array(gb, vb, v, fill));
end

function A = dimsel(A, I, dim)
% Construct the index array
szA = size(A);
szI = size(I);
stride = cumprod([1 szA(2:end)]);
if stride(dim+1) == stride(end)
    I = I(:) * stride(dim) + (1-stride(dim):0)';
elseif stride(dim) == 1
    I = I(:) + (0:stride(dim+1):stride(end)-1)';
else
    I = I(:) * stride(dim) + reshape(bsxfun(@plus, (1-stride(dim):0)', 0:stride(dim+1):stride(end)-1), [], 1);
end
% Construct the output
A = reshape(A(:,I), [szA(1) szI]);
end

function v = combine_varind(varargin)
varargin = cellfun(@var_indices, varargin, 'UniformOutput', false);
N = cellfun(@numel, varargin);
[n, m] = max(N);
v = varargin{m};
varargin = varargin(N ~= 0);
if all(N == n | N == 0) && all(cellfun(@(c) c(end) == n, varargin))
    return;
end
if all(cellfun(@(c) isequal(v, c), varargin))
    return;
end
v = unique([varargin{:}]);
end
