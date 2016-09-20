%SVM_DEMO A simple demo of SVM classification
%
%   svm_demo(N)
%
% Train different binary SVM classifiers on 2D data, and visualize the
% results.
%
%IN:
%   N - Integer number of 2D points to train and test on. Default: 500.

function svm_demo(N)
if nargin < 1
    N = 500;
end

% Create some classification training data
train_data = rand(2, N);
train_class = gt_class(train_data);

% Train 3 types of SVM
kernels = {'linear', 'intersection', 'radial'};
for a = numel(kernels):-1:1
    % Use the libsvm library
    obj(a) = libsvm();
    tic;
    % Train the SVM
    train(obj(a), train_data, train_class, struct('kernel', kernels{a}, 'probabilistic', 1, 'type', 'epsilon-svr')); % Linear
    t(1,a) = toc;
end

% Create some test data
test_data = rand(2, N);
test_class = gt_class(test_data);

% Classify test data using the 3 learned classifiers
for a = numel(obj):-1:1
    tic;
    output{a} = test(obj(a), test_data); % Linear
    t(2,a) = toc;
end

% Render the results
% Times & accuracy
for a = 1:numel(kernels)
    fprintf('Kernel: %15s. Training: %8gs. Test: %8gs. Accuracy: %g%%.\n', kernels{a}, t(1,a), t(2,a), 100*mean((test_class' < 0.5) == (output{a} < 0.5)));
end
% Classifications
figure;
subplot(221);
render(train_data, train_class);
title 'Training data'
subplot(222);
render(test_data, output{1});
title 'Linear classification'
subplot(223);
render(test_data, output{2});
title 'Intersection classification'
subplot(224);
render(test_data, output{3});
title 'Radial basis classification'
end

% Use a cubic class boundary
function class = gt_class(data)
class = ([1, -1.5, 0.57]/0.07 * [data(1,:).^3; data(1,:).^2; data(1,:)]) < data(2,:);
end

% Render the points
function render(data, prob)
hold on
patch(data(1,:), data(2,:), double(prob), 'FaceColor', 'none', 'EdgeColor', 'none', 'Marker', '+', 'MarkerEdgeColor', 'flat', 'CDataMapping', 'scaled');
hold off
caxis([0 1]);
colormap(squeeze(real2rgb(1:256, [1 0 0; 0 0 1])));
end