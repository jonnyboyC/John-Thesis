function [x_val_sort, y_val_sort, x_idx, y_idx, x_val, y_val, bnd_idx] = ...
    boundary_check_chabot(x, y, muz)
%% Find boundaries of a given image, using the mean streamwise velocity

dimensions = size(x);

x_increase = ones(dimensions);
y_increase = ones(dimensions);

% create increasing x matrix for easy logical indexing
for i = 1:dimensions(2)
    x_increase(:,i) = i*x_increase(:,i);
end

% create increasing y matrix for easy logical indexing
for i = 1:dimensions(1)
    y_increase(i,:) = i*y_increase(i,:);
end

% TODO verify that this slightly different behavior is acceptable

% Make all none zero mean velocity 1
bnd_idx = double((muz ~= 0));

% Intially make points with zero mean velocity -1
bnd_idx(bnd_idx == 0) = -1;

% Use built in edge detection to change points to 0
bnd_idx(edge(bnd_idx, 'canny')) = 0;

% Logical index points
x_val = x(bnd_idx == 0);
y_val = y(bnd_idx == 0);
x_idx = x_increase(bnd_idx == 0);
y_idx = y_increase(bnd_idx == 0);

% Put into sorted order
if size(x_val, 2) > 1
    [x_val_sort, idx] = sort(x_val);
    y_val_sort = y_val(idx);
    x_idx = x_idx(idx);
    y_idx = y_idx(idx);
else
    x_val_sort = [];
    y_val_sort = [];
    x_idx = [];
    y_idx = [];
    x_val = [];
    y_val = [];
end
