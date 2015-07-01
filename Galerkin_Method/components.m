function [pod_UdX, pod_U, vol_frac, l] = ...
    components(X, pod_U, dimensions, vol_frac, num_modes, num_elem, uniform, bnd_idx, bnd_X)
% COMPONENT, calculate the derivative terms needed to calculate the
% galerkin coefficients
% 
% type edit COMPONENT for input details

% Calculate coefficients for for pod_u's & pod_v's derivatives
[pod_UdX, pod_Ud2X] = derivatives(pod_U, bnd_idx, bnd_X, X, uniform, dimensions);

[x, u] = flow_comps_ip(X, pod_U);
dims = flow_dims(X);

% Convert all matrix quantities to vectors
for i = 1:dims;
    for j = 1:dims;
        pod_UdX.(u{i}).(x{j}) = reshape(pod_UdX.(u{i}).(x{j}), num_elem, num_modes);
        pod_Ud2X.(u{i}).(x{j}) = reshape(pod_Ud2X.(u{i}).(x{j}), num_elem, num_modes);
    end
end

% Strip data points that are superfluous
[pod_UdX, pod_Ud2X, pod_U, vol_frac] = ...
    strip_boundaries(bnd_idx, pod_UdX, pod_Ud2X, pod_U, vol_frac);

% Calculated Laplacian
for i = 1:dims
    d2pod_U.(u{i}) = 0;
    for j = 1:dims
       d2pod_U.(u{i}) = d2pod_U.(u{i}) + pod_Ud2X.(u{i}).(x{j});
    end
end

% Determine linear term
l = 0;
for i = 1:dims
    l = l + inner_prod(d2pod_U.(u{i}), pod_U.(u{i}), vol_frac);
end
end