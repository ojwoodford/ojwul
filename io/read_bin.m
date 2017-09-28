%READ_BIN Read an entire file in as an array of a given type
%
%   A = read_bin(fname, type)
%
%IN:
%   fname - string containing the filename of the file to be read.
%   type - string indicating the datatype of the values in the file.
%
%OUT:
%   A - Nx1 array of the values in the file, of datatype type.

function A = read_bin(fname, type)
fh = fopens(fname, 'r');
if fh == -1
    error('Could not open file %s for reading.', fname);
end
A = fread(fh, ['*' type]);
fclose(fh);
