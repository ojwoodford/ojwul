%WRITE_TEXT Write out an array to a text file
%
%   write_text(A, fname, append)
%
% Writes out an array to a text file using the minimum number of
% significant figures required to reconstruct the exact binary number when
% read in.
%
% The array is written out row by row, with dimensions 3 and higher
% concatenated. I.e. if the input array is 20x10x3 then the text file will
% have 20x3=60 rows, each containing 10 numbers.
%
%IN:
%   A - array to be saved.
%   fname - string containing the name of the output file.
%   append - boolean indicating whether the file is to be appended to or
%            overwritten (default).

function write_text(A, fname, append)
if ischar(A) || isstring(A)
    fmt = '%s';
else
    assert(isnumeric(A), 'Exporting non-numeric variables not supported');
    assert(isreal(A), 'Exporting complex numbers not tested');
    A = permute(A, [2 1 3]);
    A = A(:,:);
    switch class(A)
        case 'double'
            fmt = '%.16g ';
        case 'single'
            fmt = '%.8g ';
        otherwise
            fmt = '%d ';
    end
    fmt = [repmat(fmt, [1 size(A, 1)]) '\n'];
end
permission = 'wt';
if nargin > 2 && append(1)
    permission = 'at';
end
fh = fopens(fname, permission);
if fh == -1
    error('Could not open file %s for writing.', fname);
end
fprintf(fh, fmt, A);
fclose(fh);