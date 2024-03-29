function uniform = check_mesh(X)
% CHECK_MESH determine if the mesh contains uniform spacing
%
%   uniform = CHECK_MESH(X) determine if the mesh is uniform by checking
%   that the step different standard deviation is below 1e-6 that of the mean

% Get number of active dimensions
dims = flow_dims(X);

% get the field names of these active dimensions
x = flow_comps_ns(X);

% preallocate
mean_step = zeros(dims);
std_step = zeros(dims);

% create shifted matrices then take the difference and std-dev in
% difference
for i = 1:dims
    for j = 1:dims
        idx1 = flow_index({[1, -1]},i,X);
        idx2 = flow_index({[2, 0]},i,X);
        mean_step(i,j) = mean(mean(X.(x{j})(idx1{:}),i) - mean(X.(x{j})(idx2{:}),i));
        std_step(i,j) = std(mean(X.(x{j})(idx1{:}),i) - mean(X.(x{j})(idx2{:}),i));
    end
end

% if ratio between mean and std_dev is small assume mesh is uniform
if all(abs(std_step./(mean_step+eps)) < 1e-6)
    uniform = true;
else
    uniform = false;
end
end