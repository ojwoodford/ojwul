function D = cam_sqdist(P, ref)

% Compute camera centre locations
n = size(P, 3);
if size(P, 1) == 3
    for a = n:-1:1
        t(:,a) = -P(:,1:3,a) \ P(:,end,a);
    end
else
    for a = n:-1:1
        t(:,a) = (P(end,:,a) / -P(1:3,:,a))';
    end
end

if nargin > 1
    % Compute distance to reference camera
    if isscalar(ref)
        t_ = t(:,ref);
    else
        if size(P, 1) == 3
            t_ = -ref(:,1:3) \ ref(:,end);
        else
            t_ = (ref(end,:) / -ref(1:3,:))';
        end
    end
    t = bsxfun(@minus, t, t_);
    D = sum(t .* t)';
else
    % Compute distance between all cameras
    t2 = sum(t .* t);
    D = bsxfun(@plus, t2', t2) - (2 * t') * t;
end