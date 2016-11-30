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
                    case 'se3'
                        % Generators for se3
                        generators = zeros(4,4,6);
                        generators(1,4,1) = 1;
                        generators(2,4,2) = 1;
                        generators(3,4,3) = 1;
                        generators(3,2,4) = 1;
                        generators(2,3,4) = -1;
                        generators(1,3,5) = 1;
                        generators(3,1,5) = -1;
                        generators(1,2,6) = -1;
                        generators(2,1,6) = 1;
                    case 'so3'
                        % Generators for so3
                        generators = zeros(3,3,3);
                        generators(3,2,1) = 1;
                        generators(2,3,1) = -1;
                        generators(1,3,2) = 1;
                        generators(3,1,2) = -1;
                        generators(1,2,3) = -1;
                        generators(2,1,3) = 1;
                    otherwise
                        error('Lie group not recognized');
                end
            end
            this.sz = [size(generators, 1) size(generators, 2)];
            this.G = reshape(generators, [], size(generators, 3));
            this.Gv = bsxfun(@times, this.G, 1 ./ sum(abs(this.G), 1))';
        end
        
        % VEE - Convert from Lie matrix to Lie tangent space
        function tangent = vee(this, matrix)
            tangent = this.Gv * reshape(matrix, size(this.Gv, 2), []);
        end
        
        % LOG - Convert from transform to Lie matrix to Lie tangent space
        function tangent = log(this, transform)
            [~, ~, N] = size(transform);
            for a = N:-1:1
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
            for a = N:-1:1
                transform(:,:,N) = expm(this.hat(tangent(:,a)));
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
       