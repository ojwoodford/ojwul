%WRITE_BIN Write out an array to a binary file
%
%   write_bin(A, fname)
%
% Writes out an array to a binary file in the format of the data in the
% array. I.e. if the array is of type uint32 then the values are saved to
% file as 32-bit unsigned integers. This function cannot save complex
% numbers.
%
%IN:
%   A - array to be saved.
%   fname - string containing the name of the output file.

function write_bin(A, fname)
assert(isnumeric(A), 'Exporting non-numeric variables not supported');
assert(isreal(A), 'Exporting complex numbers not tested');
fh = fopen(fname, 'w', 'n');
if fh == -1
    error('Could not open file %s for writing.', fname);
end
fwrite(fh, A, class(A));
fclose(fh);