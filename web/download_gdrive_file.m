%DOWNLOAD_GDRIVE_FILE Download a file from Google Drive
%
%   download_gdrive_file(fname, fid)
%
% Download a file from Google Drive. The file must be publicly accessible.
%
%IN:
%   fname - String filename indicating where to store the download.
%   fid - String of the file's unique Google Drive ID.

% Copyright (C) Oliver Woodford 2025

function download_gdrive_file(fname, fid)
session = restful_websession('https://drive.google.com/uc');
response = session.download(fname, 'export', 'download', 'id', fid);
assert(response.StatusCode == matlab.net.http.StatusCode.OK, 'Connection failed');
% Handle the virus check confirmation (for large files)
try
    assert(isequal(response.Body.ContentType.Type, "text") && isequal(response.Body.ContentType.Subtype, "html"))
catch
    return;
end
str = read_write_entire_textfile(fname);
try
    assert(isequal('<!DOCTYPE html><html', str(1:min(end, 20))));
    assert(isequal({{'Google Drive - Virus scan warning'}}, regexp(str, '<title>(.*)</title>', 'tokens')));
    str = regexp(str, '<form id="download-form" action=".*" method="get">(.*)</form>', 'tokens');
    str = regexp(str{1}{1}, '<input type="hidden" name="([^"]*)" value="([^"]*)">', 'tokens');
    str = horzcat(str{:});
catch
    error('Response not recognized');
end
session.set_options('UseProgressMonitor', true);
delete(fname);
response = session.download(fname, str{:});
assert(response.StatusCode == matlab.net.http.StatusCode.OK, 'Download failed');
end
