function save_progress(filename, varargin)
if isempty(filename)
    return
end
% Get the global data structure
persistent save_progress_data
% Generate the unique tag name
tag = filename(isstrprop(filename, 'alphanum'));
% Check if file exists
v = exist(filename, 'file');
v = v == 0 || v == 7;
if v && isfield(save_progress_data, tag)
    % Clear the recording structure
    save_progress_data = rmfield(save_progress_data, tag);
end
% For each variable...
s = struct;
for a = 1:numel(varargin)
    % Get the iteration number
    try
        iter = save_progress_data.(tag).(varargin{a}) + 1;
    catch
        iter = 1; % First one
    end
    % Update the save iteration
    save_progress_data.(tag).(varargin{a}) = iter;
    % Capture the variable being saved and store locally in a structure
    % with the sequential name
    s.(sprintf('%s%d', varargin{a}, iter)) = evalin('caller', varargin{a});
end
if v
    % Create the file; don't compress for speed & compatibility
    save(filename, '-v6', '-struct', 's');
else
    % Append variables to the file
    save(filename, '-append', '-struct', 's');
end
return