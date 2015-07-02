function [pod_UdX, pod_U, bnd_X, volume] = ...
        components_ws(X, pod_U, dimensions, volume, num_modes, num_elem, uniform, bnd_idx, bnd_X)
% COMPONENT_WS, calculate the derivative terms needed to calculate the
% galerkin coefficients of the weak formulation
% 
% type edit COMPONENT for input details
    

% Calculate coefficients for for pod_u's & pod_v's derivatives
pod_UdX = derivatives(pod_U, bnd_idx, bnd_X, X, uniform, dimensions);

[x, u] = flow_comps_ip(X, pod_U);
dims = flow_dims(X);

for i = 1:dims;
    for j = 1:dims;
        pod_UdX.(u{i}).(x{j}) = reshape(pod_UdX.(u{i}).(x{j}), num_elem, num_modes);
    end
end

% Strip data points that are superfluous
[pod_UdX, pod_U, bnd_X, volume] = strip_boundaries(bnd_idx, pod_UdX, pod_U, bnd_X, volume);

end