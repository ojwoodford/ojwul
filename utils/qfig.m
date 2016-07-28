%QFIG  Quietly select the figure
%
%    fh = qfig(fn)
%
% Quietly selects the figure specified, without bringing it into focus
% (unless the figure doesn't exist yet).
%
% IN:
%    fn - scalar positive integer, or figure handle indicating the figure
%         to select.
%
% OUT:
%    fh - handle to the figure.

function f = qfig(f)
try
    set(0, 'CurrentFigure', f);
catch
    figure(f);
end
if nargout > 0
    f = get(0, 'CurrentFigure');
end
end