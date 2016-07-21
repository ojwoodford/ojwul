function x = whiten(x, dim)
x = bsxfun(@minus, x, mean(x, dim));
x = bsxfun(@times, x, 1./var(x, 0, dim));
end
