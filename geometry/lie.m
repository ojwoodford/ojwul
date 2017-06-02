classdef lie
    properties (Hidden = true, SetAccess = protected)
        G; % Generators for computing tangent
        Gv; % Generators for computing matrix
        sz;
    end
    methods
        function this = lie(generators)
            if ischar(generators)
                switch generators
                    case 'se2'
                        % Generators for se2
                        generators = zeros(3, 3, 3);
                        generators(1,3,1) = 1;
                        generators(2,3,2) = 1;
                        generators(1,2,3) = -1;
                        generators(2,1,3) = 1;
                    case 'aff2'
                        % Generators for aff2 (2D affine transformation)
                        generators = zeros(3, 3, 6);
                        generators(1,3,1) = 1; % X translation
                        generators(2,3,2) = 1; % Y translation
                        generators(1,2,3) = -1; % Rotation
                        generators(2,1,3) = 1;
                        generators(1,1,4) = 1; % Uniform scale
                        generators(2,2,4) = 1;
                        generators(1,1,5) = 1; % Aspect ratio
                        generators(2,2,5) = -1;
                        generators(1,2,6) = 1; % Shear
                        generators(2,1,6) = 1;
                    case 'so3'
                        generators = so3();
                    case 'rxso3'
                        generators = rxso3();
                    case 'se3'
                        % Generators for se3
                        generators = zeros(4, 4, 6);
                        generators(1,4,1) = 1;
                        generators(2,4,2) = 1;
                        generators(3,4,3) = 1;
                        generators(3,2,4) = 1;
                        generators(1:3,1:3,4:6) = so3();
                    case 'sim3'
                        % Generators for sim3
                        generators = zeros(4, 4, 7);
                        generators(1,4,1) = 1;
                        generators(2,4,2) = 1;
                        generators(3,4,3) = 1;
                        generators(3,2,4) = 1;
                        generators(1:3,1:3,4:7) = rxso3();
                    case 'sl3'
                        % Generators for sl3
                        % From "Homography-based 2D Visual Tracking and
                        % Servoing", Benhimane & Malis
                        generators = zeros(3, 3, 8);
                        generators(1,3,1) = 1; % x translation
                        generators(2,3,2) = 1; % y translation
                        generators(1,2,3) = 1;
                        generators(1:3,1:3,4:6) = so3();
                        generators(1,1,7) = 1;
                        generators(2,2,7) = -1;
                        generators(2,2,8) = -1;
                        generators(3,3,8) = 1;
                    otherwise
                        error('Lie group not recognized');
                end
            end
            this.sz = [size(generators, 1) size(generators, 2)];
            this.G = reshape(generators, [], size(generators, 3));
            this.Gv = bsxfun(@times, this.G, 1 ./ sum(abs(this.G), 1))';
        end
        
        % NDIMS - Return the dimensionality of the space
        function n = ndims(this)
            n = size(this.G, 2);
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
            transform = expm(this.hat(tangent(:,1))); % Hack for non-numeric types, e.g. autodiff
            for a = N:-1:2
                transform(:,:,a) = expm(this.hat(tangent(:,a)));
            end
        end
        
        % ADJOINT - Compute the adjoint matrix of a Lie tangent vector
        function adj = adjoint(this, tangent)
            
        end
        
        % LIEBRACKET - Apply the Lie bracket to two Lie tangent vectors
        function tangent = liebracket(this, tangentA, tangentB)
            %tangent = this.adjoint(tangentA) * tangentB;
            A = this.hat(tangentA);
            B = this.hat(tangentB);
            tangent = this.vee(A * B - B * A);
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
end

function G = so3()
% Generators for so3
G = zeros(3, 3, 3);
G(3,2,1) = 1;
G(2,3,1) = -1;
G(1,3,2) = 1;
G(3,1,2) = -1;
G(1,2,3) = -1;
G(2,1,3) = 1;
end

function G = rxso3()
% Generators for rxso3
G = so3();
G(1,1,4) = 1;
G(2,2,4) = 1;
G(3,3,4) = 1;
end
       