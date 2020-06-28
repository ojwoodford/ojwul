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
            fprintf('Value:\n');
            disp(obj.value, varargin{:});
            fprintf('Gradient:\n');
            disp(obj.deriv, varargin{:});
            fprintf('Gradient indices:\n');
            disp(obj.varind, varargin{:});
        end
        
        % Elementwise operators
        function c = plus(a, b)
            c = double(a) + double(b);
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
                ga = a.deriv .* shiftdim(db, -1);
            end
            if isautodiff(b)
                gb = b.deriv .* shiftdim(da, -1);
            end
            c = da .* db;
            [d, v] = combine_grads(ga, gb, var_indices(a), var_indices(b), size(c), 'add');
            c = autodiff(c, v, d);
        end
        
        function c = rdivide(a, b)
            da = double(a);
            db = double(b);
            ga = [];
            gb = [];
            if isautodiff(a)
                ga = a.deriv .* shiftdim(db, -1);
            end
            if isautodiff(b)
                gb = b.deriv .* -shiftdim(da, -1);
            end
            c = da ./ db;
            [d, v] = combine_grads(ga, gb, var_indices(a), var_indices(b), size(c), 'add');
            c = autodiff(c, v, d .* shiftdim(1 ./ (db .* db), -1));
        end
        
        function c = ldivide(a, b)
            c = rdivide(b, a);
        end
        
        function c = power(a, b)
            da = double(a);
            db = double(b);
            ga = [];
            gb = [];
            c = da .^ db;
            if isautodiff(a)
                ga = a.deriv .* shiftdim((da .^ (db-1)) .* db, -1);
            end
            if isautodiff(b)
                gb = b.deriv .* shiftdim(c .* log(da), -1);
            end
            [d, v] = combine_grads(ga, gb, var_indices(a), var_indices(b), size(c), 'add');
            c = autodiff(c, v, d);
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
            c = autodiff(c, a.varind, shiftdim(c, -1) .* a.deriv);
        end
        
        function c = log(a)
            c = autodiff(log(a.value), a.varind, a.deriv .* (1 ./ shiftdim(a.value, -1)));
        end
        
        function c = sqrt(a)
            c = sqrt(a.value);
            c = autodiff(c, a.varind, a.deriv .* shiftdim(min(0.5 ./ c, 1e300), -1));
        end
        
        function c = sin(a)
            c = autodiff(sin(a.value), a.varind, a.deriv .* shiftdim(cos(a.value), -1));
        end
        
        function c = cos(a)
            c = autodiff(cos(a.value), a.varind, a.deriv .* shiftdim(-sin(a.value), -1));
        end
        
        function c = tan(a)
            c = sec(a.value);
            c = autodiff(tan(a.value), a.varind, a.deriv .* shiftdim(c .* c, -1));
        end
        
        function c = asin(a)
            c = autodiff(asin(a.value), a.varind, a.deriv .* shiftdim(1 ./ sqrt(1 - a.value .* a.value), -1));
        end
        
        function c = acos(a)
            c = autodiff(acos(a.value), a.varind, a.deriv .* shiftdim(-1 ./ sqrt(1 - a.value .* a.value), -1));
        end
        
        function c = atan(a)
            c = autodiff(atan(a.value), a.varind, a.deriv .* shiftdim(1 ./ (1 + a.value .* a.value), -1));
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
            v = combine_varind(a, b);
            sz = [numel(v) size(c)];
            if isautodiff(a)
                g = reshape(reshape(grad(a, v), [], size(da, 2)) * db, sz);
            else
                g = 0;
            end
            if isautodiff(b)
                g = g + reshape(sum(shiftdim(da, -1) .* permute(grad(b, v), [1 4 2 3]), 3), sz);
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
            d = -shiftdim(c, -1) .* reshape(a.deriv, sz([1 end 2:end-1]));
            d = reshape(sum(d, 3), sz);
            d = d .* shiftdim(c, -2);
            d = reshape(sum(d, 3), sz);
            c = autodiff(c, a.varind, d);
        end
        
        function c = pinv(a)
            c = pinv(a.value);
            ct = c';
            d = permute(a.deriv, [1 3 2]);
            sz = [size(a.deriv) 1];
            if sz(2) >= sz(3)
                d_ = reshape(sum(d .* shiftdim(ct, -2), 3), sz(1), sz(3), sz(3));
                d_ = d_ - reshape(sum(shiftdim(c * a.value, -1) .* reshape(d_, sz(1), 1, sz(3), sz(3)), 3), sz(1), sz(3), sz(3));
                d_ = d_ - reshape(sum(shiftdim(c, -1) .* reshape(a.deriv, sz(1), 1, sz(2), sz(3)), 3), sz(1), sz(3), sz(3));
                d = reshape(d - reshape(sum(sum(d .* shiftdim(a.value, -2), 3) .* shiftdim(c, -3), 4), sz(1), sz(3), sz(2)), sz(1), 1, sz(3), sz(2));
                d = reshape(sum(d_ .* shiftdim(c, -2), 3), sz(1), sz(3), sz(2)) + reshape(sum(shiftdim(c * ct, -1) .* d, 3), sz(1), sz(3), sz(2));
            else
                d_ = reshape(sum(shiftdim(ct, -1) .* reshape(d, sz(1), 1, sz(3), sz(2)), 3), sz(1), sz(2), sz(2));
                d_ = d_ - reshape(sum(d_ .* shiftdim(a.value * c, -2), 3), sz(1), sz(2), sz(2));
                d_ = d_ - reshape(sum(a.deriv .* shiftdim(c, -2), 3), sz(1), sz(2), sz(2));
                d = d - reshape(sum(shiftdim(c, -1) .* reshape(sum(shiftdim(a.value, -1) .* reshape(d, sz(1), 1, sz(3), sz(2)), 3), sz(1), 1, sz(2), sz(2)), 3), sz(1), sz(3), sz(2));
                d = reshape(sum(shiftdim(c, -1) .* reshape(d_, sz(1), 1, sz(2), sz(2)), 3), sz(1), sz(3), sz(2)) + reshape(sum(d .* shiftdim(ct * c, -2), 3), sz(1), sz(3), sz(2));
            end
            c = autodiff(c, a.varind, d);
        end
        
        function c = expm(a)
            c = expm(a.value);
            sz = [size(a.deriv) 1];
            d = shiftdim(c, -1) .* reshape(a.deriv, sz([1 end 2:end-1]));
            d = reshape(sum(d, 3), sz);
            c = autodiff(c, a.varind, d);
        end
        
        function c = logm(a)
            sz = [size(a.deriv) 1];
            d = shiftdim(inv(a.value), -1) .* reshape(a.deriv, sz([1 end 2:end-1]));
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
            c = autodiff(prod(a.value, dim), a.varind, sum(a.deriv .* shiftdim(prod(c, dim) .* (1 ./ c), -1), dim+1));
        end
        
        function c = select(a, b, M)
            da = double(a);
            db = double(b);
            v = combine_varind(a, b);
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
            if k < n
                c = size(a.value, k);
            else
                sz = size(a.value);
                if numel(sz) < k
                    c = 1;
                else
                    c = prod(sz(k:end));
                end
            end
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
            v = combine_varind(a, b);
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
            v = combine_varind(varargin{:});
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
        
        function c = isfloat(a)
            c = isfloat(a.value);
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
        
        function c = isnan(a)
            c = isnan(a.value) | shiftdim(any(isnan(a.deriv)), 1);
        end
        
        function c = isfinite(a)
            c = isfinite(a.value) & shiftdim(all(isfinite(a.deriv)), 1);
        end
        
        % Other functions
        function c = ojw_interp2(I, x, y, varargin)
            if isautodiff(x) && isautodiff(y) && ~isautodiff(I)
                [c, d] = ojw_interp2(I, x.value, y.value, varargin{:});
                c = autodiff(c, x.varind, x.deriv .* d(1,:,:,:) + y.deriv .* d(2,:,:,:));
            elseif isautodiff(I) && ~isautodiff(x) && ~isautodiff(y)
                [h, w, n] = size(I.value);
                sz = size(I.value);
                [sz(1), sz(2)] = size(x);
                c = [shiftdim(I.value, -1); I.grad];
                c = reshape(reshape(c, size(c, 1), numel(I.value))', h, w, []);
                c = ojw_interp2(c, x, y, varargin{:});
                c = reshape(c, size(c, 1), size(c, 2), n, numel(I.varind)+1);
                c = autodiff(reshape(c(:,:,:,1), sz), I.varind, reshape(permute(c(:,:,:,2:end), [4 1 2 3]), [numel(I.varind) sz]));
            else
                error('Unexpected variables');
            end
        end
        
        % Alternative method avoids a subsref of autodiff variables (for
        % speed)
        function c = ojw_interp2_alt(I, x, varargin)
            if isautodiff(x)  && ~isautodiff(I)
                [c, d] = ojw_interp2(I, x.value(:,:,1), x.value(:,:,2), varargin{:});
                c = autodiff(c, x.varind, x.deriv(:,:,:,1) .* d(1,:,:,:) + x.deriv(:,:,:,2) .* d(2,:,:,:));
            else
                error('Unexpected variables');
            end
        end
        
        function c = expm_srt_3d(dP)
            dP.value(1:3,:) = dP.value(1:3,:) * sqrt(2);
            dP.deriv(:,1:3,:) = dP.deriv(:,1:3,:) * sqrt(2);
            switch size(dP, 1)
                case 3
                    type = 'so3';
                case 6
                    type = 'se3';
                    dP = autodiff(dP.value([4:6 1:3],:), dP.varind, dP.deriv(:,[4:6 1:3],:));
                case 7
                    type = 'sim3';
                    dP.value(7,:) = dP.value(7,:) * sqrt(3);
                    dP.deriv(:,7,:) = dP.deriv(:,7,:) * sqrt(3);
                    dP = autodiff(dP.value([4:6 1:3 7],:), dP.varind, dP.deriv(:,[4:6 1:3 7],:));
                otherwise
                    error('Size of dP not recognized');
            end
            c = exp(lie(type), dP);
            c = autodiff(c.value(1:3,:,:), c.varind, c.deriv(:,1:3,:,:));
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
        
        function c = proj(a)
            sz = size(a.value);
            sz(1) = sz(1) - 1;
            c = a.value(1:end-1,:);
            b = 1 ./ a.value(end,:);
            c = c .* b;
            b = shiftdim(b, -1);
            d = a.deriv(:,1:end-1,:) .* b - (shiftdim(c, -1) .* b) .* a.deriv(:,end,:);
            c = autodiff(reshape(c, sz), a.varind, reshape(d, [size(d, 1) sz]));
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
    sz2 = max(sz2, 1);
    c = repmat(gb, sz ./ sz2);
    v = vb;
    return;
end
if isempty(vb)
    sz2 = size(ga);
    sz = [sz2(1) sz];
    sz2(end+1:numel(sz)) = 1;
    sz2 = max(sz2, 1);
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
v = combine_varind_({va, vb});
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
    I = I(:) * stride(dim) + reshape((1-stride(dim):0)' + (0:stride(dim+1):stride(end)-1), [], 1);
end
% Construct the output
A = reshape(A(:,I), [szA(1) szI]);
end

function v = combine_varind(varargin)
v = combine_varind_(cellfun(@var_indices, varargin, 'UniformOutput', false));
end

function v = combine_varind_(varinds)
varinds = varinds(~cellfun(@isempty, varinds));
[n, m] = max(cellfun(@(c) c(end), varinds));
if numel(varinds{m}) == n
    return;
end
v = false(1, n);
v([varinds{:}]) = true;
v = find(v);
end
