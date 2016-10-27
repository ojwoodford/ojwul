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

% Divide by zero factor
dbzf = 1e-38;

% Default values
if nargin < 3
    s = 0;
    if nargin < 2
        t = [];
    end
end

% Angles
theta2 = r(:)' * r(:);
theta = sqrt(theta2 + dbzf);

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
    eta_r = expf(s) - (s * x + theta * y) / (s * s + theta2 + dbzf);
    eta_i = (s * y - theta * x) / (theta * (s * s + theta2) + dbzf);
end
eta_r = eta_r / (theta2 + dbzf);
A = W * eta_i + W2 * eta_r + (eye(3) * expf(s));

% Compute the translation part
M(4,4) = 1;
M(1:3,4) = A * t(:);

    function x = sinc(x)
        x = (sin(x) + dbzf) / (x + dbzf);
    end

    function x = cosf(x)
        x = ((1 - cos(x)) + dbzf) / (x * x + dbzf);
    end

    function x = expf(x)
        x = ((exp(x) - 1) + dbzf) / (x + dbzf);
    end
end
