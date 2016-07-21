function h = render_coordinate_frame(P)
%RENDER_COORDINATE_FRAME Renders RGB coordinate frames given by P

% Invert the frames
P(4,4,:) = 1;
for a = 1:size(P, 3)
    P(:,:,a) = inv(P(:,:,a));
end
P = P(1:3,:,:);

% Plot the axes coordinates
colour = {[0.8 0 0], [0 0.6 0], [0 0 0.9]};
rshp = @(x) permute(x, [2 3 1]);
O = rshp(P(:,4,:));
N = NaN(size(O));
hold on
for a = 3:-1:1
    X = reshape(cat(1, O, O+rshp(P(:,a,:)), N), [], 3);
    h(a) = plot3(X(:,1), X(:,2), X(:,3), 'r-', 'Color', colour{a});
end
end