%LOAD_FIELD Load only one field from a file
%
%   x = load_field(name, field)
%
%IN:
%   name - Filename string of the file containing the field.
%   field - Name string of the field to be loaded.
%OUT:
%   x - Loaded field.

function x = load_field(name, field)
x = load(name, field);
x = x.(field);
end