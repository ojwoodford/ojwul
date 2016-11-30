function saveh5(filename, dataAsStruct)
%SAVEH5 Saves a struct dataset into a HDF5 format file
%
% Examples:
%    str.a=1;str.b=2;
%    saveh5(filename, str);
%
% This function saves data to a HDF5 file. Each member of the root struct
% is saved a seperate dataset in the file.
%
% IN:
%   filename - A string containing the path to a HDF5 file.
%   dataAsStruct - A structure with data in members.
% NOTE: no nested support

fields = fieldnames(dataAsStruct);
numFields = numel(fields);
%Just in case
if exist(filename, 'file')
    delete(filename);
end

for f = 1:numFields
    %Drop trailing 1s
    sz = size(dataAsStruct.(fields{f}));
    last1 = find(sz ~= 1,1,'last');
    if isempty(last1)
        sz = 1;
    else
        sz = sz(1:find(sz ~= 1,1,'last'));
    end
    h5create(filename, ['/' fields{f}], sz, 'Datatype', class(dataAsStruct.(fields{f})));
    h5write(filename, ['/' fields{f}], dataAsStruct.(fields{f}));
end
