%READ_FLOAT32 Read an entire file in as an array of 32-bit floats
%
%   A = read_float32(fname)
%
%IN:
%   fname - string containing the filename of the file to be read.
%
%OUT:
%   A - Nx1 single array of the values in the file.

function A = read_float32(fname)
fh = fopen(fname, 'r');
A = fread(fh, '*float32');
fclose(fh);