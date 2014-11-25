function [pod_2D] = regroup(pod, dimensions)
%% Transform 2D matrix of ensembles back to sets of 2D images
pod_modes = size(pod,2);
x_dim = dimensions(1);
y_dim = dimensions(2);

pod_2D = zeros(x_dim, y_dim, pod_modes);

% Convert from 2D to 3D
for i = 1:pod_modes
    pod_2D(1:x_dim,1:y_dim,i) = reshape(pod(:,i), dimensions);
end

