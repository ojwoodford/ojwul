%EXPM_SRT_3D Compute a transformation matrix, given the Lie algebra vector
%
%   M = expm_srt_3d_sym(r, [t, [s]])
%
% Computes the symbolic transformation matrix defined by a Lie vector
% consisting of a rotation, translation and uniform scaling.
%
% This function applies the formula given in the paper:
% "Distances and Means of Direct Similarities"
% M-T Pham et al.
%
%IN:
%   r - 3x1 rotation parameter vector
%   t - 3x1 translation parameter vector. Default: [].
%   s - scalar uniform scaling parameter. Default: 0.
%
%OUT:
%   M - 3x3 (if t == []) or 4x4 transformation matrix.

function M = expm_srt_3d_sym(r, t, s)

% Default values
if nargin < 3
    s = sym(0);
    if nargin < 2
        t = [];
    end
end

% Angles
theta2 = r(:)' * r(:);
theta = sqrt(theta2);

% Skew matrix
W = skew(r);
W2 = W * W;

% Rotation matrix via Rodrigues formula
M = W * sinc(theta) + W2 * cosf(theta) + eye(3);

% Apply scale
M = M * exp(s);

if isempty(t)
    return;
end

% Translation multipler
if nargin < 3
    eta_r = sym(1) - sinc(theta);
    eta_i = cosf(theta);
else
    x = cos(theta) * exp(s) - 1;
    y = sin(theta) * exp(s);
    eta_r = expf(s) - (s * x + theta * y) / (s * s);
    eta_i = (s * y - theta * x) / (theta * (s * s + theta2));
end
eta_r = at0is1(theta2, @(theta2) eta_r / theta2);
A = W * eta_i + W2 * eta_r + (eye(3) * expf(s));

% Compute the translation part
M(1:3,4) = simplify(A * t(:));
end

function x = sinc(x)
x = at0is1(x, @(x) sin(x) / x);
end

function x = cosf(x)
x = at0is1(x, @(x) (1 - cos(x)) / (x * x));
end

function x = expf(x)
x = at0is1(x, @(x) (exp(x) - 1) / x);
end

function x = at0is1(x, fun)
if isequal(x, sym(0))
    x = sym(1);
else
    x = fun(x);
end
end
