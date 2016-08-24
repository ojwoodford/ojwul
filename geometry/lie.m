classdef lie
    properties (Hidden = true, SetAccess = protected)
        G; % Generators
        normalizers;
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
                    otherwise
                        error('Lie group not recognized');
                end
            end
            this.sz = [size(generators, 1) size(generators, 2)];
            this.G = reshape(generators, [], size(generators, 3));
            this.normalizers = 1 ./ sum(abs(this.G), 1)';
        end
        function tangent = log(this, transform)
            tangent = sum(bsxfun(@times, reshape(logm(transform), [], 1), this.G), 1)' .* this.normalizers;
        end
        function transfom = exp(this, tangent)
            transfom = expm(reshape(sum(bsxfun(@times, this.G, tangent'), 2), this.sz));
        end
    end
end
       