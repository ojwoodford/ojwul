%MOV2IM Convert a movie file to a sequence of images 
%
% Examples:
%   mov2im infile outfile_format
%
%IN:
%   infile - string containing the name of the input video.
%   outfile_format - format string for the movie frames. The filename for
%                    frame N is given by sprintf(outfile_format, N).

% Copyright (C) Oliver Woodford 2015

function mov2im(infile, outfile_format, varargin)

% Create the movie object
hMov = VideoReader(infile);

% Write out the frames
frame = 0;
while hasFrame(hMov)
    imwrite(readFrame(hMov), sprintf(outfile_format, frame), varargin{:});
    frame = frame + 1;
end
end
