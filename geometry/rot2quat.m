%ROT2QUAT Convert an array of 3x3 rotation matrices into quaternions
%
%    Q = rot2quat(R, [w_last])
%
% This function generates quaternions of the form [w x y z]', such that R =
% eye(3) corresponds to Q = [1 0 0 0]'. If w_last is true, quaternions are
% of the form [x y z w].
%
%IN:
%   R - 3x3xN array of rotation matrices.
%   w_last - boolean indicating whether quaternions should be of the form
%            [x y z w]'. Default: false.
%
%OUT:
%   Q - 4xN array of corresponding quaternions.

function Q = rot2quat(R, w_last)
M = cat(3, [1 0 0 0 0 0 0; 0 0 0 0 1 0 -1; 0 0 -1 0 0 1 0; 0 1 0 -1 0 0 0], ...
           [0 0 0 0 1 0 -1; 1 0 0 0 0 0 0; 0 1 0 1 0 0 0; 0 0 1 0 0 1 0], ...
           [0 0 -1 0 0 1 0; 0 1 0 1 0 0 0; 1 0 0 0 0 0 0; 0 0 0 0 1 0 1], ...
           [0 1 0 -1 0 0 0; 0 0 1 0 0 1 0; 0 0 0 0 1 0 1; 1 0 0 0 0 0 0]);
sz = size(R);
R = reshape(R, 9, []);
[m, n] = max([1, 1, 1; 1, -1, -1; -1, 1, -1; -1, -1, 1] * R([1 5 9],:));
m = sqrt(1 + m) / 2;
Q = sum(bsxfun(@times, M(:,:,n), shiftdim([m; bsxfun(@times, R([2 3 4 6 7 8],:), 1 ./ (4 * m))], -1)), 2);
Q = reshape(Q, [4 sz(3:end) 1]);
if nargin > 1 && w_last
    Q = Q([2:4 1],:);
end
end
