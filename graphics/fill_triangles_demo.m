%FILL_TRIANGLES_DEMO  Demonstrates fill_triangles
%
%   fill_triangles_demo([resolution, [use_single]])
%
% Demonstrates usage of the fill_triangles function by rendering a
% cube with flat shaded, gouraud shaded and texture mapped sides.
%
% Rotate the cube using arrow keys. Quit by pressing q.
%
% IN:
%     resolution - resolution of scene to render, in pixels. Default: 200.
%     use_single - non-zero if single precision is to be used, otherwise
%                  (default) double precision is used. Single precision
%                  approximately doubles the speed of the triangle
%                  renderer.

function fill_triangles_demo(resolution, use_single)

clf reset;
%set(gcf,'DoubleBuffer','on');
map = jet(256);
map = map(end:-1:1,:);
colormap(map);
map = [];

if nargin < 1
    resolution = 400;
end

% Cube vertices
CP = [-1 -1 -1 ; 1 -1 -1 ; 1 1 -1 ; -1 1 -1; ...
      -1 -1  1 ; 1 -1  1 ; 1 1  1 ; -1 1  1]';

% Cube vertex colours for gouraud shaded faces
CC = colorcube(8)' * 255;

% Cube triangles (must be clockwise outwards)
T = [1 3 2 ; 1 4 3 ; ...
     1 2 5 ; 2 6 5 ; ...
     2 3 6 ; 3 7 6 ; ...
     3 4 7 ; 4 8 7 ; ...
     1 8 4 ; 1 5 8; ...
     5 6 7 ; 5 7 8]';

% Initialise parameters
I = zeros(4, resolution, resolution);
I(1,:,:) = 1/6;
Q = zeros(6, 36);
half_res = resolution / 2;
perspective_correct = 1;

% Initialise the gouraud shaded faces
P = CC(:,T);

% Initialise the flat shaded faces
P(:,7:12) = repmat([255 0 0]', [1 6]); % Red
P(:,19:24) = repmat([0 0 255]', [1 6]); % Blue

% Initialise the texture mapped faces
% Face 1
TM1 = double(imread('cameraman.tif'));
TM1 = TM1(:,:,[1 1 1]);
P(1,1:6) = 350; % Index of the texture
TV = [1 1; size(TM1, 2) 1; size(TM1, 2) size(TM1, 1); 1 size(TM1, 1)]';
P(2:3,1:6) = TV(:,T(1:6)); % Coordinates on the texture
% Face 1
TM2 = double(imread('peppers.png'));
P(1,31:36) = 450; % Index of the texture
TV = [1 1 ; size(TM2, 2) 1 ; size(TM2, 2) size(TM2, 1); 1 size(TM2, 1)]';
P(2:3,31:36) = TV(:,T(31:36)-4); % Coordinates on the texture
P = repmat(P, [2 1]);
P = reshape(P, [6 3 12]);

one = 1;
if nargin > 1 && use_single
    I = single(I);
    P = single(P);
    TM1 = single(TM1);
    TM2 = single(TM2);
    one = single(1);
end

time = zeros(1, 6);
f = 0;

%R = eye(3);
R = [cos(pi/3) -sin(pi/3) 0 ; sin(pi/3) cos(pi/3) 0; 0 0 1] * ...
    [cos(-pi/6) 0 -sin(-pi/6) ; 0 1 0 ; sin(-pi/6) 0 cos(-pi/6)] * ...
    [1 0 0 ; 0 cos(pi/6) -sin(pi/6) ; 0 sin(pi/6) cos(pi/6)];
incr = pi / 60; % 3 degrees
Rx = [1 0 0 ; 0 cos(incr) -sin(incr) ; 0 sin(incr) cos(incr)];
Ry = [cos(incr) 0 -sin(incr) ; 0 1 0 ; sin(incr) 0 cos(incr)];
lims = [5 - sqrt(3) 6];
lims(2) = 255 / (lims(2) - lims(1));

% Set up the button callback for interactive control
fig_h = gcf;
figure(fig_h);
pos = get(0, 'ScreenSize');
pos = round(pos/2);
pos = [pos(3:4) resolution*2+300 resolution+300];
pos(1:2) = pos(1:2) - round(pos(3:4)/2);
set(fig_h, 'KeyPressFcn', @keypress_callback, 'Units', 'pixels', 'Position', pos);

% Create the plot layout
subplot(121);
H1 = imshow(zeros(resolution, resolution, 3, 'uint8'));
set(gca, 'xlimmode', 'manual', 'ylimmode', 'manual',...
         'zlimmode', 'manual', 'climmode', 'manual',...
         'alimmode', 'manual');
title({'Cube with 2 flat shaded', '2 gouraud shaded & 2', 'texture mapped sides'});
xlabel({'Use arrow keys', 'to rotate the cube.', 'Press q to quit.', 'Press p to toggle', 'perspective correct', 'interpolation'});

subplot(122);
H2 = imshow(zeros(resolution, resolution, 3, 'uint8'));
set(gca, 'xlimmode', 'manual', 'ylimmode', 'manual',...
         'zlimmode', 'manual', 'climmode', 'manual',...
         'alimmode', 'manual');
title('Calculated depth map of cube.')

while true
    t = tic;
    
    X = R * CP;
    X(3,:) = X(3,:) + 5;
    X(1,:) = resolution * X(1,:) ./ X(3,:) + half_res;
    X(2,:) = resolution * X(2,:) ./ X(3,:) + half_res;
    
    time(6) = time(6) + toc(t);
    t = tic;
    
    % Cull back facing triangles
    Q = reshape(X(1:2,T), [2 3 12]);
    Q(:,1,:) = Q(:,2,:) - Q(:,1,:);
    Q(:,2,:) = Q(:,3,:) - Q(:,2,:);
    Q = squeeze((Q(1,1,:) .* Q(2,2,:)) - (Q(2,1,:) .* Q(1,2,:)));
    L = find(Q < 0);
    Q = reshape(P(:,:,L), [6 3*length(L)]);
    L = reshape(T(:,L), [1 3*length(L)]);
    Q(1:3,:) = X(:,L);
        
    time(5) = time(5) + toc(t);
    t = tic;
    
    if perspective_correct
        % Ensure perspective correct interpolation by dividing through by depth
        V = 1 ./ Q(3,:);
        Q(3,:) = V;
        for a = 4:size(Q, 1)
            Q(a,:) = Q(a,:) .* V;
        end
        
        time(2) = time(2) + toc(t);
        t = tic;
        
        % Call the triangle renderer
        J = fill_triangles(I, Q, -1);

        time(1) = time(1) + toc(t);
        t = tic;
        
        % Divide through by interpolated inverse depth, for perspective correct
        % interpolation
        V = one ./ J(1,:);
        J(1,:) = V;
        for a = 2:size(J, 1)
            J(a,:) = J(a,:) .* V;
        end
        time(2) = time(2) + toc(t);
    else
        % Call the triangle renderer
        J = fill_triangles(I, Q, 0);
        time(1) = time(1) + toc(t);
    end
    t = tic;
    
    % Render the texture mapped surfaces
    L = J(2,:,:) > 400;
    if any(L(:))
        J(2:4,L) = squeeze(ojw_interp2(TM2, J(3,L), J(4,L), 'l', one))';
    else
        L = J(2,:,:) > 300 & J(2,:,:) <= 400;
        if any(L(:))
            J(2:4,L) = squeeze(ojw_interp2(TM1, J(3,L), J(4,L), 'l', one))';
        end
    end
    time(3) = time(3) + toc(t);
    t = tic;

    % Display the output
    J = shiftdim(J, 1);
    Z = J(:,:,1);
    RGB = uint8(J(:,:,2:end));
    
    set(H1, 'CData', RGB);
    Z = uint8((Z - lims(1)) * lims(2)) + 1;
    set(H2, 'Cdata', Z);
    drawnow;
    time(4) = time(4) + toc(t);
    f = f + 1;
    
    % Get the key press
    while true
        set(fig_h, 'userData', '')
        waitfor(fig_h, 'userData')
        button = lower(get(fig_h, 'userData'));
        switch button
            case 31
                R = Rx * R;
                break
            case 30
                R = Rx' * R;
                break
            case 29
                R = Ry * R;
                break
            case 28
                R = Ry' * R;
                break
            case 'p'
                perspective_correct = xor(1, perspective_correct);
                if perspective_correct
                    I(1,:,:) = 1/6;
                else
                    I(1,:,:) = 6;
                end
                break
            case {'q', '', 13}
                close(gcf);
                total = sum(time);
                fprintf('Rendering speed of %0.0f frames per second, with processing time allocated as follows:\n', f / total);
                total = total - time(4);
                time = time * 100;
                fprintf('   Projection -             %0.1f%%\n', time(6) / total);
                fprintf('   Back face culling -      %0.1f%%\n', time(5) / total);
                fprintf('   Triangle rendering -     %0.1f%%\n', time(1) / total);
                fprintf('   Perspective correction - %0.1f%%\n', time(2) / total);
                fprintf('   Texture mapping -        %0.1f%%\n', time(3) / total);
                fprintf('And drawing to the screen accounting for %0.1f times all these time put together!\n', time(4) / (100 * total));
                return
        end
    end
end
return

function keypress_callback(source_handle, event_data)
set(source_handle, 'userData', event_data.Character);
return
    

