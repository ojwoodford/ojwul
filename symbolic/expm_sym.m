%EXPM_SYM Symbolic matrix exponential, up to a certain order
%
%   B = expm_sym(A, order)
%
% Given a symbolic (or numeric) matrix, this function computes an
% approximation of the matrix exponential up to a certain order.
%
%IN:
%   A - MxM input matrix.
%   order - scalar indicating the order up to which the matrix exponential
%           is to be computed.
%
%OUT:
%   B - MxM output matrix.

function B = expm_sym(A, order)
B = eye(size(A));
if isa(A, 'sym')
    B = sym(B);
end
A_ = B;
for a = 1:order
    A_ = A_ * (A / a);
    B = B + A_;
end