%ROT2QUAT Convert an array of 3x3 rotation matrices into quaternions
%
% Q = rot2quat(R)

function Q = rot2quat(R)

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
