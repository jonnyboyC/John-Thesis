function [pod_2D] = regroup(pod, dimensions)
pod_modes = size(pod,2);
x_dim = dimensions(1);
y_dim = dimensions(2);

pod_2D = zeros(x_dim, y_dim, pod_modes);

for i = 1:pod_modes
    pod_2D(1:x_dim,1:y_dim,i) = reshape(pod(:,i), dimensions);
end

