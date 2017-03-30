%RNG_SEEDER Seed the random number generator, and print the seed if generated
%
%   seed = rng_seeder()
%   rng_seeder(seed)
%
% This function intializes the random number generator, and prints out the
% seed if one is not given or output.
%
%IN:
%   seed - scalar seed for the random number generator.
%
%OUT:
%   seed - scalar seed used to initialize the random number generator.

function seed = rng_seeder(seed)
if nargin < 1 || isempty(seed)
    rng('shuffle');
    seed = ceil(rand(1) * (2^31));
    if nargout < 1
        fprintf('RNG seed: %d\n', seed);
    end
end
rng(seed);
end
