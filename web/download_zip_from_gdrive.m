function download_zip_from_gdrive(name, fid)
tname = tempname();
co = onCleanup(@() delete_dir(tname));
mkdir(tname);
download_unzip(tname, fid, name);
end

function delete_dir(tname)
if exist(tname, 'dir')
    rmdir(tname, 's');
end
end

function download_unzip(tname, fid, name)
base = cd();
temp_cd(tname);
fname = strcat(name, '.zip');
session = restful_websession('https://drive.google.com/uc');
response = session.download(fname, 'export', 'download', 'id', fid);
assert(response.StatusCode == matlab.net.http.StatusCode.OK, 'Connection failed');
try
    unzip(fname);
catch
    % Handle no virus check confirmation (for large files)
    str = read_write_entire_textfile(fname);
    delete(fname);
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
    response = session.download(fname, str{:});
    assert(response.StatusCode == matlab.net.http.StatusCode.OK, 'Download failed');
    unzip(fname);
end
delete(fname);
if exist(name, 'dir')
    movefile(name, base);
else
    cd(base);
    movefile(tname, name);
end
end
