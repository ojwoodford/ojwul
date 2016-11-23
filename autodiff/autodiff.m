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
                    obj.deriv = reshape(eye(numel(a)), [numel(a) size(a)]);
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
            end
            c = zeros([numel(vars) size(obj.value)]);
            [Lia, Locb] = ismember(vars, obj.varind);
            c(Lia,:) = obj.deriv(Locb,:);
        end
        
        function c = var_indices(obj)
            c = obj.varind;
        end
        
        % Elementwise operators
        function c = plus(a, b)
            [c, v] = combine_grads(grad(a), grad(b), var_indices(a), var_indices(b), 'add');
            c = autodiff(bsxfun(@plus, double(a), double(b)), v, c);
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
            [c, v] = combine_grads(ga, gb, var_indices(a), var_indices(b), 'add');
            c = autodiff(bsxfun(@times, da, db), v, c);
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
            [c, v] = combine_grads(ga, gb, var_indices(a), var_indices(b), 'add');
            c = autodiff(bsxfun(@rdivide, da, db), v, bsxfun(@times, c, shiftdim(1 ./ (db .* db), -1)));
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
            da = double(a);
            db = double(b);
            if isautodiff(b)
                v = b.varind;
                c = sum(bsxfun(@times, shiftdim(da, -1), permute(b.deriv, [1 4 2 3])), 3);
                if isautodiff(a)
                    c = c + sum(bsxfun(@times, a.deriv, shiftdim(db, -2)), 3);
                end
                c = reshape(c, [size(b.deriv, 1) size(da, 1) size(db, 2)]);
            else
                v = a.varind;
                c = sum(bsxfun(@times, a.deriv, shiftdim(db, -2)), 3);
                c = reshape(c, [size(a.deriv, 1) size(da, 1) size(db, 2)]);
            end
            c = autodiff(da * db, v, c);
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
        
        % Reduction methods: Sum, prod, min, max
        function c = sum(a, dim)
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
            if isscalar(da)
                da = repmat(da, size(db));
            end
            if isscalar(db)
                da(M) = db;
            else
                da(M) = db(M);
            end
            if isautodiff(a)
                ga = reshape(a.deriv, [], numel(da));
                if isscalar(da)
                    ga = repmat(ga, [1 numel(db)]);
                end
                if isautodiff(b)
                    gb = reshape(b.deriv, [], numel(db));
                    if isscalar(db)
                        ga(:,M) = repmat(gb, [1, sum(M(:))]); 
                    else
                        ga(:,M) = gb(:,M);
                    end
                else
                    ga(:,M) = 0;
                end
            else
                ga = reshape(b.deriv, [], numel(db));
                if isscalar(db)
                    ga = repmat(ga, [1 numel(da)]);
                end
                ga(:,~M) = 0;
            end
            c = autodiff(da, v, reshape(ga, [size(ga, 1) size(da)]));
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
                c = autodiff(c, v, dimsel(a.deriv, d, dim));
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
                c = autodiff(c, v, dimsel(a.deriv, d, dim));
            end
        end
        
        % Access functions
        function c = end(a, k, n)
            c = size(a.value, k);
        end
         
        function c = subsref(a, s)
            assert(strcmp(s.type, '()'));
            s_ = [{':'} s.subs];
            c = autodiff(a.value(s.subs{:}), a.varind, a.deriv(s_{:}));
        end
        
        function c = subsasgn(a, s, b)
            assert(strcmp(s.type, '()'));
            c = double(a);
            c(s.subs{:}) = double(b);
            d = grad(a);
            s_ = [{':'} s.subs];
            d(s_{:}) = grad(b);
            c = autodiff(c, a.varind, d);
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
                v_ = {[size(a.deriv, 1) varargin{1}]};
            else
                v_ = [{size(a.deriv, 1)}; varargin(:)];
            end
            c = autodiff(reshape(a.value, varargin{:}), a.varind, reshape(a.deriv, v_{:}));
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
        
        % Other functions        
        function c = ojw_interp2(I, x, y, interp_mode, oobv)
            if nargin < 5
                oobv = NaN;
            end
            assert(nargin < 4 || interp_mode(1) == 'l', 'Only linear interpolation supported');
            assert(isautodiff(x) && isautodiff(y) && ~isautodiff(I), 'Unexpected variables');
            [c, d] = ojw_interp2(I, x.value, y.value, 'l', oobv);
            c = autodiff(c, x.varind, bsxfun(@times, x.deriv, d(1,:,:,:)) + bsxfun(@times, y.deriv, d(2,:,:,:)));
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

function [c, v] = combine_grads(ga, gb, va, vb, M)
if ischar(M)
    switch M
        case 'add'
            func = @plus;
            fill = 0;
        case 'mult'
            func = @times;
            fill = 1;
    end
    M = [];
else
    func = @plus;
    fill = 0;
end
if isempty(va)
    c = gb;
    v = vb;
    return;
end
if isempty(vb)
    c = ga;
    v = va;
    return;
end
if isequal(va, vb)
    c = bsxfun(func, ga, gb);
    v = va;
    return;
end
v = unique(va, vb);
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
