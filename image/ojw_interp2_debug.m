%OJW_INTERP2_DEBUG Image sampling, with plotting of the sample points
function varargout = ojw_interp2_debug(varargin)
% Forward the call, for error checking
[varargout{1:nargout}] = ojw_interp2(varargin{:});
% Render the sample points
figure(7429);
clf reset;
imdisp(varargin{1});
hold on;
plot(double(varargin{2}), double(varargin{3}), 'r.');
end