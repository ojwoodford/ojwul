function h = ojw_texture_map(Z, I)
%OJW_TEXTURE_MAP  Plot a texture mapped surface
%
%   h = ojw_texture_map(Z, I)
%
% Texture maps a regular mesh in a Matlab figure. If Z and I are the same
% size then the texture map is linerly interpolated, otherwise it is mapped
% as is.
%
%IN:
%   Z - MxN[x3] Matrix of depths, or regular grid of 3d points.
%   I - Texture map image.
%
%OUT:
%   h - Handle to the surface object

% $Id: ojw_texture_map.m,v 1.3 2008/11/20 13:36:08 ojw Exp $

[m n c] = size(Z);
if c == 3
    X = Z(:,:,1);
    Y = Z(:,:,2);
    Z = Z(:,:,3);
else
    X = 1:n;
    Y = 1:m;
    Z = Z(:,:,1);
end
params = {'edgecolor', 'none'};
if nargin > 1
    if numel(I) < 4
        h = surface(X, Y, Z, 'facecolor', I, params{:});
        light;
        lighting phong;
    else
        [h w c] = size(I);
        mx = max(reshape(I, h*w*c, 1));
        I = double(I);
        if mx > 1
            I = I / 255;
        end
        if isequal([h w], [m n])
            h = surface(X, Y, Z, I, 'facecolor', 'interp', params{:});
        else
            h = surface(X, Y, Z, I, 'facecolor', 'texturemap', params{:});
        end
    end
else
    h = surface(X, Y, Z, 'facecolor', 'interp', params{:});
end
if nargout < 1
    clear h
end
return