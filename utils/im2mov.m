%IM2MOV Convert a sequence of images to a movie file
%
% Examples:
%   im2mov infile outfile
%   im2mov(A, outfile)
%   im2mov(..., '-fps', n)
%   im2mov(..., '-quality', q)
%   im2mov(..., '-profile', profile) 
%   im2mov(..., '-nocrop')
%   
% This function converts an image sequence to a movie.
%
% To create a video from a series of figures, export to an image sequence
% using export_fig, then convert to a movie, as follows:
%
%    frame_num = 0;
%    for a = 2 .^ (3:6)
%       peaks(a);
%       export_fig(sprintf('test.%3.3d.png', frame_num), '-nocrop');
%       frame_num = frame_num + 1;
%    end
%    im2mov('test.000.png', 'output', '-fps', 0.5, '-profile', 'MPEG-4');
%
%IN:
%   infile - string containing the name of the first input image.
%   outfile - string containing the name of the output video (without an
%             extension).
%   A - HxWxCxN array of input images, stacked along fourth dimension, to
%       be converted to a movie.
%   -nocrop - option indicating that the borders of the output are not to
%             be cropped.
%   -fps - option pair, the value of which gives the number of frames per
%          second. Default: 15.
%   -quality - option pair, the value of which gives the video quality,
%              from 0 to 100. Default: 75.
%   -profile - option pair, the value of which is a profile string passed
%              to VideoWriter object. Default: 'MPEG-4'.

% Copyright (C) Oliver Woodford 2013

function im2mov(A, varargin)

% Parse the input arguments
[A, options] = parse_args(A, varargin{:});

% Get a consistent interface to the sequence
if isnumeric(A)
    im = @(a) A(:,:,:,a);
    num_frames = size(A, 4);
elseif ischar(A)
    im = imstream(A, 50);
    num_frames = im.num_frames;
end

crop_func = @(A) A;
if options.crop ~= 0
    %% Determine the borders
    A = im(1);
    [h, w, c] = size(A);
    L = w;
    R = w;
    T = h;
    B = h;
    bcol = A(1,1,:);
    for n = num_frames:-1:1
        A = im(n);
        [h, w, c] = size(A);
        bail = false;
        for l = 1:min(L, w)
            for a = 1:c
                if ~all(col(A(:,l,a)) == bcol(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
        L = l;
        bail = false;
        for r = min(w, R):-1:L
            for a = 1:c
                if ~all(col(A(:,r,a)) == bcol(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
        R = r;
        bail = false;
        for t = 1:min(T, h)
            for a = 1:c
                if ~all(col(A(t,:,a)) == bcol(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
        T = t;
        bail = false;
        for b = min(h, B):-1:T
            for a = 1:c
                if ~all(col(A(b,:,a)) == bcol(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
        B = b;
    end
    crop_func = @(A) A(T:B,L:R,:);
end

% Create the movie object
hMov = VideoWriter(options.outfile, options.profile);
set(hMov, 'FrameRate', options.fps);
try
    % Not all profiles support quality
    set(hMov, 'Quality', options.quality);
catch
end
open(hMov);
% Write out the frames
for a = 1:num_frames
    writeVideo(hMov, crop_func(im(a)));
end
close(hMov);
return

%% Parse the input arguments
function [infile, options] = parse_args(infile, varargin)
% Set the defaults
options = struct('outfile', '', ...
                 'crop', true, ....
                 'profile', 'MPEG-4', ...
                 'quality', 75, ...
                 'fps', 15);

% Go through the arguments
a = 0;
n = numel(varargin);
while a < n
    a = a + 1;
    if ischar(varargin{a}) && ~isempty(varargin{a})
        if varargin{a}(1) == '-'
            opt = lower(varargin{a}(2:end));
            switch opt
                case 'nocrop'
                    options.crop = false;
                otherwise
                    if ~isfield(options, opt)
                        error('Option %s not recognized', varargin{a});
                    end
                    a = a + 1;
                    if ischar(varargin{a}) && ~ischar(options.(opt))
                        options.(opt) = str2double(varargin{a});
                    else
                        options.(opt) = varargin{a};
                    end
            end
        else
            options.outfile = varargin{a};
        end
    end
end

if isempty(options.outfile)
    if ~ischar(infile)
        error('No output filename given.');
    end
    % Generate the output filename from the input filename
    [path, outfile] = fileparts(infile);
    options.outfile = fullfile(path, outfile);
end
return

function A = col(A)
A = A(:);
return
