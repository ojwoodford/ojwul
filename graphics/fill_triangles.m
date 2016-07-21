%FILL_TRIANGLES  Standard, linear interpolation triangle renderer
%
%   fill_triangles(I, P, zbuff)
%   J = fill_triangles(I, P, zbuff)
%
% A software implementation of a standard, linear interpolation (i.e.
% incorrect for perspective), z-buffered triangle renderer.
%
% IN:
%     I - CxMxN double or single input image array. A single array makes the
%         rendering faster, but less accurate. C is the number of channels
%         to be interpolated. If there are only two dimensions then C is
%         assumed to be 1, and M and N the given dimensions. If zbuff != 0
%         then the first channel must be the initialised z-buffer. The
%         remaining channels are the initial image colours or texture
%         vertices, depending on what is being interpolated.
%     P - (2+C)x(3T) matrix, of the same type as I. In the first two
%         rows are the x and y image coordinates respectively, and the rest
%         of the rows are the channel values of each vertex, in the same
%         order as in I. The number of columns must be a multiple of 3,
%         each group of 3 columns being the 3 verices of a triangle, with T
%         triangles in total. Note that shared vertices have multiple
%         entries, allowing for flat shaded faces.
%     zbuff - 1, 0 or -1. If 0, z-buffering is not used, and pixels are
%             simply painted over each other in the order of triangles
%             given. 1 means that smaller z-depths are deemed nearer, and
%             -1 means that smaller z-depths are deemed further away.
%
% OUT:
%     J - CxMxN output image array, of same type as I. Each channel is the
%         same as that from I, but updated according to the list of
%         triangles provided, and the z-test for each pixel (if z-buffering
%         is used). WARNING: If no output argument is provided, the input
%         array I is written into directly; use with care!

function varargout = fill_triangles(varargin)
sourceList = {'fill_triangles.cpp'}; % Cell array of source files
[varargout{1:nargout}] = compile(varargin{:}); % Compilation happens here
return