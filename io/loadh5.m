function data = loadh5(filename)
%LOADH5 Loads data from a HDF5 format file
%
% Examples:
%    X = loadh5(filename)
%
% This function loads data from a HDF5 file and returns each item as a 
% member of a structure.
% Note: no support for nested items
%
% IN:
%   filename - A string containing the path to a HDF5 file.
%
% OUT:
%   data - A structure containing each database in the HDF5 file as a
%          member.

if ~exist(filename, 'file')
    error('The specified file does not exist:\n%s', filename);
end

%Load info
fileInfo = hdf5info(filename);

data = struct();

for i = 1:numel(fileInfo.GroupHierarchy.Datasets)
    buffer = hdf5read(filename, fileInfo.GroupHierarchy.Datasets(i).Name);
    
    %Handle arrays
    if isa(buffer,'hdf5.h5string')
        if numel(buffer) > 1
            data = setfield(data, fileInfo.GroupHierarchy.Datasets(i).Name(2:end), arrayfun(@(x) x.data, buffer,'UniformOutput',false)); %#ok<SFLD>
        else 
            data = setfield(data, fileInfo.GroupHierarchy.Datasets(i).Name(2:end), buffer.data);
        end
    elseif isa(buffer,'hdf5.h5array')
        data = setfield(data, fileInfo.GroupHierarchy.Datasets(i).Name(2:end), arrayfun(@(x) x.data, buffer,'UniformOutput',false)); %#ok<SFLD>
    else
        data = setfield(data, fileInfo.GroupHierarchy.Datasets(i).Name(2:end), buffer); %#ok<SFLD>
    end
end
