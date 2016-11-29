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
                    case 'so3'
                        generators = so3();
                    case 'se3'
                        % Generators for se3
                        generators = zeros(4, 4, 6);
                        generators(1,4,1) = 1;
                        generators(2,4,2) = 1;
                        generators(3,4,3) = 1;
                        generators(3,2,4) = 1;
                        generators(1:3,1:3,4:6) = so3();
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
            if N == 1
                % Hack for non-numeric types, e.g. autodiff
                transform = expm(this.hat(tangent));
                return;
            end
            for a = N:-1:1
                transform(:,:,N) = expm(this.hat(tangent(:,a)));
            end
        end
        
        % LIEBRACKET - Apply the Lie bracket to two Lie tangent vectors
        function tangent = liebracket(this, tangentA, tangentB)
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
       