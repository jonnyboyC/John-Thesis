function [bnd_x, bnd_y, bnd_idx] = boundary_check_chabot(x, mean_u)
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
bnd_idx = double((mean_u ~= 0));

% Intially make points with zero mean velocity -1
bnd_idx(bnd_idx == 0) = -1;

% Use built in edge detection to change points to 0
bnd_idx(edge(bnd_idx, 'canny')) = 0;

% Determine the vector normal of the boundary
[bnd_x, bnd_y] = gradient(bnd_idx);
