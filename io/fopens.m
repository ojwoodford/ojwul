%FOPENS  Open file which is always closed when the function exits.
%   FID = FOPENS(FILENAME) opens the file FILENAME for read access.
%   FILENAME is a string containing the name of the file to be opened.  
%   (On PC systems, FOPENS opens files for binary read access.)
%
%   FILENAME can be a MATLABPATH relative partial pathname.  If the
%   file is not found in the current working directory, FOPENS searches for 
%   it on the MATLAB search path.  On UNIX systems, FILENAME may also start
%   with a "~/" or a "~username/", which FOPENS expands to the current
%   user's home directory or the specified user's home directory,
%   respectively.
%
%   FID is a scalar MATLAB integer valued double, called a file identifier. 
%   You use FID as the first argument to other file input/output
%   routines, such as FREAD and FCLOSE. If FOPENS cannot open the file, it 
%   throws an error.
%
%   FID = FOPENS(FILENAME,PERMISSION) opens the file FILENAME in the
%   mode specified by PERMISSION:
%   
%       'r'     open file for reading
%       'w'     open file for writing; discard existing contents
%       'a'     open or create file for writing; append data to end of file
%       'r+'    open (do not create) file for reading and writing
%       'w+'    open or create file for reading and writing; discard 
%               existing contents
%       'a+'    open or create file for reading and writing; append data 
%               to end of file
%       'W'     open file for writing without automatic flushing
%       'A'     open file for appending without automatic flushing
%   
%   FILENAME can be a MATLABPATH relative partial pathname only if the file
%   is opened for reading.
%
%   You can open files in binary mode (the default) or in text mode.
%   In binary mode, no characters get singled out for special treatment.
%   In text mode on the PC, the carriage return character preceding
%   a newline character is deleted on input and added before the newline
%   character on output.  To open a file in text mode, append 't' to the
%   permission string, for example 'rt' and 'w+t'.  (On Unix, text and
%   binary mode are the same, so this has no effect.  On PC systems
%   this is critical.)
%
%   If the file is opened in update mode ('+'), you must use an FSEEK or 
%   FREWIND between an input command like FREAD, FSCANF, FGETS, or FGETL 
%   and an output command like FWRITE or FPRINTF.  You must also use an 
%   FSEEK or FREWIND between an output command and an input command.
%
%   Two file identifiers are automatically available and need not be
%   opened.  They are FID=1 (standard output) and FID=2 (standard error).
%   
%   [FID, MESSAGE] = FOPENS(FILENAME,...) returns a system dependent error 
%   message if the open is not successful.
%
%   [FID, MESSAGE] = FOPENS(FILENAME,PERMISSION,MACHINEFORMAT) opens the
%   specified file with the specified PERMISSION and treats data read
%   using FREAD or data written using FWRITE as having a format given
%   by MACHINEFORMAT. MACHINEFORMAT is one of the following strings:
%
%   'native'      or 'n' - local machine format - the default
%   'ieee-le'     or 'l' - IEEE floating point with little-endian
%                          byte ordering
%   'ieee-be'     or 'b' - IEEE floating point with big-endian
%                          byte ordering
%   'ieee-le.l64' or 'a' - IEEE floating point with little-endian
%                          byte ordering and 64 bit long data type
%   'ieee-be.l64' or 's' - IEEE floating point with big-endian byte
%                          ordering and 64 bit long data type.
%   
%   [FID, MESSAGE] = FOPENS(FILENAME,PERMISSION,MACHINEFORMAT,ENCODING) 
%   opens the specified file using the specified PERMISSION and 
%   MACHINEFORMAT. ENCODING is a string that specifies the character 
%   encoding scheme associated with the file. It must be the empty 
%   string ('') or a name or alias for an encoding scheme. Some examples 
%   are 'UTF-8', 'latin1', 'US-ASCII', and 'Shift_JIS'. For common names 
%   and aliases, see the Web site 
%   http://www.iana.org/assignments/character-sets. If ENCODING is 
%   unspecified or is the empty string (''), MATLAB's default encoding 
%   scheme is used.
%
%   [FILENAME,PERMISSION,MACHINEFORMAT,ENCODING] = FOPENS(FID) returns the 
%   filename, permission, machine format, and character encoding values
%   used by MATLAB when it opened the file associated with identifier FID.
%   MATLAB does not determine these output values by reading information 
%   from the opened file. For any of these parameters that were not 
%   specified when the file was opened, MATLAB returns its default value. 
%   The ENCODING string is a standard character encoding scheme name that 
%   may not be the same as the ENCODING argument used in the call to FOPEN 
%   that opened the file. An invalid FID returns empty strings for all 
%   output arguments.
%
%   FIDS = FOPENS('all') returns a row vector containing the file 
%   identifiers for all the files currently opened by the user 
%   (but not 1 or 2).
%   
%   The 'W' and 'A' permissions do not automatically perform a flush 
%   of the current output buffer after output operations.
%
%   Any file opened with FOPENS will automatically be closed when the
%   function exits, even through an error, or when the workspace is
%   cleared.
%   
%   See also FOPEN.

function varargout = fopens(varargin)
% Silly check
if nargout == 0
    return
end

% Call fopen
[varargout{1:nargout}] = fopen(varargin{:});

% Check if opened
if isequal(varargout{1}, -1)
    % Not opened. Throw an error.
    error('File %s could not be opened', varargin{1});
else
    % Set the cleanup function
    [name, name] = fileparts(tempname);
    name = ['CloseFileOnExitObj_' name];
    assignin('caller', name, onCleanup(@() fclose_quiet(varargout{1})));
end
end

function fclose_quiet(fid)
try
    fclose(fid);
catch
end
end
