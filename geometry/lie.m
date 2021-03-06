classdef lie
    properties (Hidden = true, SetAccess = protected)
        G; % Generators for computing tangent
        Gv; % Generators for computing matrix
        sz;
    end
    methods
        function this = lie(generators_)
            if ischar(generators_)
                generators_ = generators(generators_);
            end
            this.sz = [size(generators_, 1) size(generators_, 2)];
            this.G = reshape(generators_, [], size(generators_, 3));
            this.Gv = pinv(this.G);
        end
        
        % NDIMS - Return the dimensionality of the space
        function n = ndims(this)
            n = size(this.G, 2);
        end
        
        % SIZE - Return the size of the warp
        function sz = size(this)
            sz = this.sz;
        end
        
        % VEE - Convert from Lie matrix to Lie tangent space
        function tangent = vee(this, matrix)
            tangent = this.Gv * reshape(matrix, size(this.Gv, 2), []);
        end
        
        % LOG - Convert from transform to Lie matrix to Lie tangent space
        function tangent = log(this, transform)
            [~, ~, N] = size(transform);
            tangent = this.vee(logm(transform(:,:,1))); % Hack for non-numeric types, e.g. autodiff
            for a = N:-1:2
                tangent(:,a) = this.vee(logm(transform(:,:,a)));
            end
        end
        
        % HAT - Convert from Lie tangent space to Lie matrix
        function matrix = hat(this, tangent)
            matrix = reshape(this.G * tangent, this.sz(1), this.sz(2), []);
        end
        
        % EXP - Convert from Lie tangent space to transform
        function transform = exp(this, tangent)
            [~, N] = size(tangent);
            if N == 1
                transform = expm(this.hat(tangent));
                return;
            end
            transform = expm(this.hat(tangent(:,1))); % Hack for non-numeric types, e.g. autodiff
            for a = N:-1:2
                transform(:,:,a) = expm(this.hat(tangent(:,a)));
            end
        end
        
        % ADJOINT - Compute the adjoint matrix of a Lie tangent vector
        function adj = adjoint(this, tangent)
            
        end
        
        function tangent = Adjoint_times(this, transform, tangent)
            tangent = tmult(transform, this.hat(tangent));
            for a = 1:size(transform, 3)
                transform(:,:,a) = inv(transform(:,:,a));
            end
            tangent = this.vee(tmult(tangent, transform));
        end
        
        % LIEBRACKET - Apply the Lie bracket to two Lie tangent vectors
        function tangent = liebracket(this, tangentA, tangentB)
            %tangent = this.adjoint(tangentA) * tangentB;
            A = this.hat(tangentA);
            B = this.hat(tangentB);
            tangent = this.vee(tmult(A, B) - tmult(B, A));
        end
        
        % INTERP1 - Linear interpolation between transforms
        function Vq = interp1(this, X, V, Xq)
            for a = numel(X)-1:-1:1
                tangents(:,a) = this.log(V(:,:,a+1) / V(:,:,a));
            end
            [~, ind] = histc(Xq, [X(:); Inf]);
            ind = min(max(ind, 1), numel(X)-1);
            for a = numel(Xq):-1:1
                Vq(:,:,a) = this.exp(tangents(:,ind(a)) * ((Xq(a) - X(ind(a))) / (X(ind(a)+1) - X(ind(a))))) * V(:,:,ind(a));
            end
        end
    end
    methods (Static)
        function check_orthogonal(groups)
            if nargin < 1
                groups = {'so2', 'se2', 'sim2', 'aff2', 'so3', 'rxso3', 'uv3', 'rxuv3', 'se3', 'sim3', 'sl3'};
            elseif ~iscell(groups)
                groups = {groups};
            end
            for group = groups
                % Chech the generators are orthogonal
                G = generators(group{1});
                G = reshape(G, [], size(G, 3));
                dot_prod = tril(squeeze(sum(bsxfun(@times, G, permute(G, [1 3 2])), 1)), -1);
                [y, x] = find(abs(dot_prod) > 1e-14);
                if ~isempty(y)
                    y = [y(:)'; x(:)'];
                    warning('The following pairs of %s group generators are not orthogonal:\n%s', group{1}, sprintf('   %d,%d\n', y));
                end
            end
        end
    end
end

function G = generators(group)
switch group
    case 'shift2'
        % Generators for 2D translations
        G = zeros(3, 3, 2);
        G(1,3,1) = 1; % X translation
        G(2,3,2) = 1; % Y translation
    case 'so2'
        % Generators for so2
        G = [0 -1; 1 0];
    case 'se2'
        % Generators for se2
        G = generators('shift2');
        G(1:2,1:2,3) = generators('so2'); % Rotation
    case 'sim2'
        % Generators for sim2
        G = cat(3, generators('se2'), diag([0.5 0.5 -1])); % Uniform scale
    case 'rs2'
        G = generators('sim2');
        G = G(:,:,3:4);
    case 'aff2'
        % Generators for aff2 (2D affine transformation)
        G = zeros(3, 3, 6);
        G(:,:,1:4) = generators('sim2');
        G(1,1,5) = 1; % Aspect ratio
        G(2,2,5) = -1;
        G(1,2,6) = 1; % Shear
        G(2,1,6) = 1;
    case 'so3'
        % Generators for so3
        G = zeros(3, 3, 3);
        G(3,2,1) = 1;
        G(2,3,1) = -1;
        G(1,3,2) = 1;
        G(3,1,2) = -1;
        G(1,2,3) = -1;
        G(2,1,3) = 1;
    case 'rxso3'
        % Generators for rxso3
        G = generators('so3');
        G(1,1,4) = 1;
        G(2,2,4) = 1;
        G(3,3,4) = 1;
    case 'uv3' % 3D unit vector
        G = generators('so3');
        G = G(:,:,2:3);
    case 'rxuv3' % 3D vector parameterized by rotation and scale
        G = generators('rxso3');
        G = G(:,:,2:4);
    case 'se3'
        % Generators for se3
        G = zeros(4, 4, 6);
        G(1,4,1) = 1;
        G(2,4,2) = 1;
        G(3,4,3) = 1;
        G(1:3,1:3,4:6) = generators('so3');
    case 'sim3'
        % Generators for sim3
        G = zeros(4, 4, 7);
        G(1,4,1) = 1;
        G(2,4,2) = 1;
        G(3,4,3) = 1;
        G(1:3,1:3,4:7) = generators('rxso3');
    case 'sl3'
        % Generators for sl3
        G = zeros(3, 3, 8);
        G_ = generators('so3');
        G(:,:,1:3) = G_;
        G_ = abs(G_);
        G(:,:,4) = 0.5 * G_(:,:,3) - G_(:,:,1) - G_(:,:,2);
        G(:,:,5) = G_(:,:,2) + G_(:,:,3) - 0.5 * G_(:,:,1);
        G(:,:,6) = G_(:,:,3) + G_(:,:,1) - 0.5 * G_(:,:,2);
        G(1,1,7) = -1;
        G(3,3,7) = 1;
        G(1,1,8) = 0.5;
        G(3,3,8) = 0.5;
        G(2,2,8) = -1;
    otherwise
        error('Lie group %s not recognized', group);
end
% Normalize the generators
G = bsxfun(@times, G, 1 ./ sqrt(sum(sum(G .* G, 1), 2)));
end
       