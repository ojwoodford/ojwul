%BITCOUNT Count the number of set bits in each column of the input
%
%   B = bitcount(A)
%
% Count the number of set bits in each column of the input array,
% typecast as a bit vector.
%
%IN:
%   A - MxNx... input array.
%
%OUT:
%   B - 1xNx... output array of bit counts.

function A = bitcount(A)

persistent lufun lutable
if isempty(lutable)
    % Generate the lookup table
    lutable = uint8(sum(dec2bin(0:255) - '0', 2));
    % Use intlut function if possible
    lufun = which('intlut');
    if isempty(lufun)
        lufun = @intlutsub;
    else
        cwd = cd();
        try
            cd([fileparts(lufun) '/private']);
            intlutmex(uint8(0), uint8(0));
            lufun = str2func('intlutmex');
        catch
            lufun = @intlut;
        end
        cd(cwd);
    end
end

% Convert to an index into the lookup table
sz = size(A);
sz(1) = 1;
[n, n] = size(A);
A = reshape(typecast(A(:), 'uint8'), [], n);

% Look up the number of set bits for each byte
A = lufun(A, lutable);

% Sum the number of set bits per column
if size(A, 1) < 32
    A = sum(A, 1, 'native');
else
    A = sum(A, 1);
end
A = reshape(A, sz);
end

function A = intlutsub(A, lutable)
A = lutable(uint16(A) + uint16(1));
end
