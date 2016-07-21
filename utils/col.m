function x = col(x, dim)
x = x(:);
if nargin > 1
    x = shiftdim(x, 1-dim);
end
end