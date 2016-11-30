%GRAD Dummy helper function for autodiff
function c = grad(a, vars)
if nargin < 2
    c = 0;
else
    c = zeros([numel(vars) size(a)]);
end
end
